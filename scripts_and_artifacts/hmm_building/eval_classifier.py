#!/usr/bin/env python3
"""
Evaluate tubulin classifier on raw UniProt sequences (not used in HMM training).
This is a proper holdout evaluation since raw sequences were filtered during balancing.
"""

import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from loguru import logger

from lib.hmm import get_raw_fasta_path
from lib.hmm.classifier import TubulinClassifier
from lib.types import TubulinFamily


@dataclass
class EvalResult:
    seq_id: str
    true_family: TubulinFamily
    predicted_family: TubulinFamily | None
    correct: bool
    top_score: float
    second_score: float | None
    score_margin: float | None


def load_raw_sequences() -> list[tuple[str, str, TubulinFamily]]:
    """Load all sequences from raw FASTAs."""
    all_seqs = []
    family_counts = {}
    
    for family in TubulinFamily:
        fasta_path = get_raw_fasta_path(family.value)
        if not fasta_path.exists():
            logger.warning(f"Missing: {fasta_path}")
            continue
        
        records = list(SeqIO.parse(fasta_path, "fasta"))
        family_counts[family] = len(records)
        
        for rec in records:
            all_seqs.append((rec.id, str(rec.seq), family))
    
    print("Raw sequences per family:")
    for family, count in family_counts.items():
        print(f"  {family.value}: {count}")
    print(f"  Total: {len(all_seqs)}")
    
    # Shuffle for random processing order
    random.seed(42)
    random.shuffle(all_seqs)
    
    return all_seqs


def evaluate(
    test_seqs: list[tuple[str, str, TubulinFamily]],
    classifier: TubulinClassifier
) -> list[EvalResult]:
    """Run classifier on all test sequences."""
    results = []
    n = len(test_seqs)
    
    for i, (seq_id, sequence, true_family) in enumerate(test_seqs):
        result = classifier.classify(
            sequence,
            rcsb_id=f"RAW_{true_family.value}",
            auth_asym_id=seq_id[:15],
            force=True
        )
        
        predicted = result.assigned_family
        correct = predicted == true_family
        
        sorted_hits = sorted(result.hits, key=lambda h: h.score, reverse=True)
        top_score = sorted_hits[0].score if sorted_hits else 0
        second_score = sorted_hits[1].score if len(sorted_hits) > 1 else None
        margin = (top_score - second_score) if second_score else None
        
        results.append(EvalResult(
            seq_id=seq_id,
            true_family=true_family,
            predicted_family=predicted,
            correct=correct,
            top_score=top_score,
            second_score=second_score,
            score_margin=margin
        ))
        
        # Progress every 100
        if (i + 1) % 100 == 0:
            current_acc = sum(1 for r in results if r.correct) / len(results)
            print(f"  {i + 1}/{n} ({current_acc:.1%} accuracy so far)")
    
    return results


def print_report(results: list[EvalResult]):
    """Print detailed evaluation report."""
    
    total = len(results)
    correct = sum(1 for r in results if r.correct)
    accuracy = correct / total if total > 0 else 0
    
    print("\n" + "=" * 80)
    print("RAW SEQUENCE EVALUATION REPORT")
    print("=" * 80)
    
    print(f"\nOverall: {correct}/{total} correct ({accuracy:.2%})")
    
    # Per-family
    print("\nPer-family breakdown:")
    print("-" * 80)
    print(f"{'Family':12s} {'Correct':>8s} {'Total':>8s} {'Accuracy':>10s} {'Avg Margin':>12s} {'Min Margin':>12s}")
    print("-" * 80)
    
    for family in TubulinFamily:
        fam_results = [r for r in results if r.true_family == family]
        if not fam_results:
            continue
        
        fam_correct = sum(1 for r in fam_results if r.correct)
        fam_acc = fam_correct / len(fam_results)
        margins = [r.score_margin for r in fam_results if r.score_margin is not None]
        avg_margin = sum(margins) / len(margins) if margins else 0
        min_margin = min(margins) if margins else 0
        
        print(f"{family.value:12s} {fam_correct:8d} {len(fam_results):8d} {fam_acc:10.2%} {avg_margin:12.1f} {min_margin:12.1f}")
    
    # Confusion matrix
    print("\nConfusion Matrix:")
    print("-" * 80)
    
    families = list(TubulinFamily)
    confusion = defaultdict(lambda: defaultdict(int))
    
    for r in results:
        pred = r.predicted_family.value if r.predicted_family else "NONE"
        confusion[r.true_family.value][pred] += 1
    
    # Header

    header_label = "True \\ Pred"
    print(f"{header_label:12s}", end="")
    for f in families:
        print(f"{f.value[:6]:>8s}", end="")
    print(f"{'NONE':>8s}")
    
    for true_fam in families:
        print(f"{true_fam.value:12s}", end="")
        for pred_fam in families:
            count = confusion[true_fam.value].get(pred_fam.value, 0)
            print(f"{count:8d}", end="")
        none_count = confusion[true_fam.value].get("NONE", 0)
        print(f"{none_count:8d}")
    
    # Misclassifications
    misclassified = [r for r in results if not r.correct]
    
    if misclassified:
        print(f"\nMisclassifications ({len(misclassified)} total):")
        print("-" * 80)
        
        # Group by error type
        error_types = defaultdict(list)
        for r in misclassified:
            pred = r.predicted_family.value if r.predicted_family else "NONE"
            error_types[(r.true_family.value, pred)].append(r)
        
        for (true, pred), errors in sorted(error_types.items(), key=lambda x: -len(x[1])):
            print(f"\n  {true} -> {pred}: {len(errors)} cases")
            # Show first few examples
            for r in errors[:3]:
                margin_val = r.score_margin if r.score_margin else 0
                print(f"    {r.seq_id[:30]:30s}  score={r.top_score:.1f}  margin={margin_val:.1f}")
            if len(errors) > 3:
                print(f"    ... and {len(errors) - 3} more")
    
    # Score statistics
    print("\nScore margin distribution:")
    print("-" * 80)
    
    correct_margins = [r.score_margin for r in results if r.correct and r.score_margin]
    incorrect_margins = [r.score_margin for r in results if not r.correct and r.score_margin]
    
    if correct_margins:
        correct_margins.sort()
        p5 = correct_margins[len(correct_margins) // 20]
        p50 = correct_margins[len(correct_margins) // 2]
        print(f"  Correct:   n={len(correct_margins):4d}  min={min(correct_margins):7.1f}  p5={p5:7.1f}  median={p50:7.1f}  max={max(correct_margins):7.1f}")
    
    if incorrect_margins:
        incorrect_margins.sort()
        p50 = incorrect_margins[len(incorrect_margins) // 2] if incorrect_margins else 0
        print(f"  Incorrect: n={len(incorrect_margins):4d}  min={min(incorrect_margins):7.1f}  median={p50:7.1f}  max={max(incorrect_margins):7.1f}")
    
    # Suggest threshold
    if incorrect_margins and correct_margins:
        # Find margin where we'd filter out most errors
        sorted_incorrect = sorted(incorrect_margins, reverse=True)
        if sorted_incorrect:
            suggested = sorted_incorrect[0] + 10  # slightly above max incorrect margin
            filtered_correct = sum(1 for m in correct_margins if m >= suggested)
            print(f"\n  Suggested margin threshold: {suggested:.0f}")
            print(f"    Would keep {filtered_correct}/{len(correct_margins)} correct ({filtered_correct/len(correct_margins):.1%})")
            print(f"    Would filter all {len(incorrect_margins)} misclassifications")
    
    print("\n" + "=" * 80)


def main():
    print("Loading raw sequences...")
    test_seqs = load_raw_sequences()
    
    print("\nInitializing classifier...")
    classifier = TubulinClassifier(use_cache=False)
    
    print("\nRunning evaluation (this may take a minute)...")
    results = evaluate(test_seqs, classifier)
    
    print_report(results)


if __name__ == "__main__":
    main()