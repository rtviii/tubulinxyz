#!/usr/bin/env python3
"""
Deep dive into misclassified sequences to understand why.
"""

from Bio import SeqIO

from lib.hmm import get_raw_fasta_path
from lib.hmm.tubulin_classifier import TubulinClassifier
from lib.models.types_tubulin import TubulinFamily


def analyze_delta_misclassifications():
    """Analyze the delta -> beta/alpha misclassifications."""
    
    classifier = TubulinClassifier(use_cache=False)
    
    # Load delta sequences
    delta_path = get_raw_fasta_path("delta")
    delta_seqs = list(SeqIO.parse(delta_path, "fasta"))
    
    print(f"Analyzing {len(delta_seqs)} delta sequences...\n")
    
    misclassified = []
    low_margin = []
    
    for rec in delta_seqs:
        result = classifier.classify(
            str(rec.seq),
            rcsb_id="DELTA",
            auth_asym_id=rec.id[:15],
            force=True
        )
        
        if result.assigned_family != TubulinFamily.DELTA:
            misclassified.append((rec, result))
        else:
            # Check margin
            sorted_hits = sorted(result.hits, key=lambda h: h.score, reverse=True)
            if len(sorted_hits) > 1:
                margin = sorted_hits[0].score - sorted_hits[1].score
                if margin < 50:  # low confidence
                    low_margin.append((rec, result, margin))
    
    # Report misclassified
    print("=" * 80)
    print(f"MISCLASSIFIED DELTA SEQUENCES ({len(misclassified)})")
    print("=" * 80)
    
    for rec, result in misclassified:
        print(f"\n{'─' * 80}")
        print(f"Sequence: {rec.id}")
        print(f"Length: {len(rec.seq)}")
        print(f"Assigned: {result.assigned_family.value if result.assigned_family else 'NONE'}")
        print()
        print(result.report())
        print()
        print(result.domain_report())
        print()
        
        # Show sequence characteristics
        seq = str(rec.seq)
        print(f"Sequence preview: {seq[:60]}...")
        print(f"                  ...{seq[-60:]}")
    
    # Report low-margin correct classifications
    print("\n" + "=" * 80)
    print(f"LOW-MARGIN DELTA CLASSIFICATIONS ({len(low_margin)})")
    print("=" * 80)
    
    for rec, result, margin in sorted(low_margin, key=lambda x: x[2])[:10]:
        print(f"\n{rec.id}: margin={margin:.1f}")
        sorted_hits = sorted(result.hits, key=lambda h: h.score, reverse=True)
        for hit in sorted_hits[:3]:
            print(f"  {hit.family.value}: {hit.score:.1f}")


def compare_hmm_profiles():
    """Show HMM model sizes to understand why confusion happens."""
    classifier = TubulinClassifier()
    
    print("\nHMM Profile Comparison:")
    print("-" * 40)
    
    info = classifier.hmm_info()
    for family, data in info.items():
        print(f"  {family:10s}: M={data['M']:4d} (from {data['nseq']} seqs)")
    
    print("\nNote: M = model length (match states)")
    print("Delta M=615 vs Beta M=448 - delta is longer")
    print("Misclassifications may be truncated/partial sequences")


if __name__ == "__main__":
    compare_hmm_profiles()
    print()
    analyze_delta_misclassifications()