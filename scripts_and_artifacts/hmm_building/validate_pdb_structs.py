#!/usr/bin/env python3
"""
Validate HMM classifier against PDB structure entities.
Compares HMM classification to existing family assignments and
analyzes non-tubulin proteins (MAPs, MIPs, etc.)

Usage:
    python validate_pdb_structures.py run          # Run analysis and save report
    python validate_pdb_structures.py show         # Show existing report
    python validate_pdb_structures.py show --file other_report.json
"""

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from datetime import datetime

from rich.console import Console
from rich.table import Table
from rich.panel import Panel
from rich import print as rprint

from lib.hmm.tubulin_classifier import TubulinClassifier, ClassificationResult
from lib.models.types_tubulin import TubulinFamily

# Paths
TUBETL_DATA = Path("./TUBETL_DATA")
DEFAULT_REPORT_PATH = Path("./PDB_Classification_Report.json")

console = Console()

@dataclass
class EntityResult:
    rcsb_id: str
    entity_id: str
    description: str
    sequence_length: int
    existing_family: str | None
    predicted_family: str | None
    confidence: str
    best_score: float
    margin: float
    all_scores: dict[str, float]
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, d: dict) -> "EntityResult":
        return cls(**d)


def get_confidence(result: ClassificationResult) -> tuple[str, float]:
    """Determine confidence level and margin."""
    if not result.hits:
        return "no_hit", 0.0
    
    sorted_hits = sorted(result.hits, key=lambda h: h.score, reverse=True)
    best = sorted_hits[0]
    margin = (best.score - sorted_hits[1].score) if len(sorted_hits) > 1 else best.score
    
    # Check coverage
    if best.domains:
        total_covered = sum(d.env_to - d.env_from + 1 for d in best.domains)
        coverage = total_covered / result.sequence_length
    else:
        coverage = 1.0
    
    if best.score < 100:
        return "not_tubulin", margin
    elif coverage < 0.7:
        return "partial", margin
    elif margin < 50:
        return "low_margin", margin
    elif best.score < 250:
        return "weak", margin
    else:
        return "high", margin


def categorize_non_tubulin(description: str) -> str:
    """Categorize non-tubulin protein by description."""
    desc_lower = (description or "").lower()
    
    if "map" in desc_lower or "microtubule-associated" in desc_lower:
        return "MAP"
    elif "mip" in desc_lower or "inner protein" in desc_lower:
        return "MIP"
    elif "stathmin" in desc_lower:
        return "Stathmin"
    elif "kinesin" in desc_lower:
        return "Kinesin"
    elif "dynein" in desc_lower:
        return "Dynein"
    elif "tau" in desc_lower:
        return "Tau"
    elif "eb1" in desc_lower or "eb3" in desc_lower or "mapre" in desc_lower:
        return "EB proteins"
    elif "ttl" in desc_lower or "ligase" in desc_lower:
        return "TTL/Ligases"
    elif "chaperone" in desc_lower or "cct" in desc_lower or "tric" in desc_lower:
        return "Chaperones"
    elif "katanin" in desc_lower or "spastin" in desc_lower or "severing" in desc_lower:
        return "Severing enzymes"
    elif "colchicine" in desc_lower or "drug" in desc_lower:
        return "Drug-related"
    elif "gtp" in desc_lower or "gtpase" in desc_lower:
        return "GTPase-related"
    elif "actin" in desc_lower:
        return "Actin"
    elif "antibody" in desc_lower or "nanobody" in desc_lower or "fab" in desc_lower:
        return "Antibody/Nanobody"
    elif "darpin" in desc_lower:
        return "DARPin"
    else:
        return "Other"


def load_structure(json_path: Path) -> dict | None:
    """Load structure JSON."""
    try:
        with open(json_path) as f:
            return json.load(f)
    except Exception as e:
        console.print(f"[red]Error loading {json_path}: {e}[/red]")
        return None


def extract_polypeptide_entities(structure: dict) -> list[dict]:
    """Extract polypeptide entities from structure."""
    entities = []
    for entity_id, entity in structure.get("entities", {}).items():
        if entity.get("type") == "polymer" and entity.get("polymer_type") == "Protein":
            entities.append({
                "entity_id": entity_id,
                "description": entity.get("pdbx_description", "Unknown"),
                "sequence": entity.get("one_letter_code_can", ""),
                "family": entity.get("family"),
                "sequence_length": entity.get("sequence_length", 0),
            })
    return entities


def run_analysis() -> dict:
    """Run full analysis and return report dict."""
    console.print("[bold blue]Running PDB Structure Validation[/bold blue]\n")
    
    classifier = TubulinClassifier(use_cache=False)
    
    all_results: list[EntityResult] = []
    
    structure_dirs = sorted([d for d in TUBETL_DATA.iterdir() if d.is_dir()])
    console.print(f"Found {len(structure_dirs)} structures\n")
    
    for i, struct_dir in enumerate(structure_dirs):
        rcsb_id = struct_dir.name
        json_path = struct_dir / f"{rcsb_id}.json"
        
        if not json_path.exists():
            continue
        
        structure = load_structure(json_path)
        if not structure:
            continue
        
        entities = extract_polypeptide_entities(structure)
        
        for entity in entities:
            sequence = entity["sequence"]
            if not sequence or len(sequence) < 20:
                continue
            
            result = classifier.classify(
                sequence,
                rcsb_id=rcsb_id,
                auth_asym_id=entity["entity_id"],
                force=True
            )
            
            confidence, margin = get_confidence(result)
            predicted = result.assigned_family.value if result.assigned_family else None
            existing = entity["family"]
            
            all_scores = {h.family.value: h.score for h in result.hits}
            best_score = max(all_scores.values()) if all_scores else 0
            
            entity_result = EntityResult(
                rcsb_id=rcsb_id,
                entity_id=entity["entity_id"],
                description=entity["description"] or "Unknown",
                sequence_length=len(sequence),
                existing_family=existing,
                predicted_family=predicted,
                confidence=confidence,
                best_score=best_score,
                margin=margin,
                all_scores=all_scores,
            )
            all_results.append(entity_result)
        
        if (i + 1) % 50 == 0:
            console.print(f"  Processed {i + 1}/{len(structure_dirs)} structures...")
    
    # Build report
    report = {
        "generated_at": datetime.now().isoformat(),
        "total_structures": len(structure_dirs),
        "total_entities": len(all_results),
        "results": [r.to_dict() for r in all_results],
    }
    
    return report


def compute_stats(results: list[EntityResult]) -> dict:
    """Compute statistics from results."""
    stats = {
        "total": len(results),
        "by_confidence": defaultdict(int),
        "by_predicted_family": defaultdict(int),
        "by_existing_family": defaultdict(int),
        "matches": 0,
        "mismatches": 0,
        "no_existing": 0,
        "non_tubulin_categories": defaultdict(int),
    }
    
    for r in results:
        stats["by_confidence"][r.confidence] += 1
        
        if r.predicted_family:
            stats["by_predicted_family"][r.predicted_family] += 1
        
        if r.existing_family:
            stats["by_existing_family"][r.existing_family] += 1
        
        if r.confidence in ("not_tubulin", "no_hit"):
            cat = categorize_non_tubulin(r.description)
            stats["non_tubulin_categories"][cat] += 1
        
        if r.existing_family and r.predicted_family:
            if r.existing_family == r.predicted_family:
                stats["matches"] += 1
            else:
                stats["mismatches"] += 1
        elif r.predicted_family and not r.existing_family:
            stats["no_existing"] += 1
    
    return stats


def display_report(report: dict):
    """Display report with rich formatting."""
    results = [EntityResult.from_dict(r) for r in report["results"]]
    stats = compute_stats(results)
    
    # === HEADER ===
    console.print(Panel(
        f"[bold]PDB Classification Report[/bold]\n"
        f"Generated: {report['generated_at']}\n"
        f"Structures: {report['total_structures']} | Entities: {report['total_entities']}",
        style="blue"
    ))
    
    # === SUMMARY STATS ===
    console.print("\n[bold cyan]═══ SUMMARY STATISTICS ═══[/bold cyan]\n")
    
    console.print("[bold]By confidence level:[/bold]")
    confidence_order = ["high", "weak", "low_margin", "partial", "not_tubulin", "no_hit"]
    for conf in confidence_order:
        count = stats["by_confidence"].get(conf, 0)
        pct = count / stats["total"] * 100 if stats["total"] > 0 else 0
        color = {
            "high": "green",
            "weak": "yellow",
            "low_margin": "yellow",
            "partial": "magenta",
            "not_tubulin": "cyan",
            "no_hit": "red",
        }.get(conf, "white")
        bar = "█" * int(pct / 2)
        console.print(f"  [{color}]{conf:15s}[/{color}]: {count:4d} ({pct:5.1f}%) {bar}")
    
    console.print("\n[bold]Predicted tubulin families:[/bold]")
    for family in ["alpha", "beta", "gamma", "delta", "epsilon"]:
        count = stats["by_predicted_family"].get(family, 0)
        console.print(f"  {family:10s}: {count:4d}")
    
    console.print("\n[bold]Existing vs Predicted:[/bold]")
    console.print(f"  [green]Matches[/green]:     {stats['matches']}")
    console.print(f"  [red]Mismatches[/red]:  {stats['mismatches']}")
    console.print(f"  [yellow]No existing[/yellow]: {stats['no_existing']}")
    
    # === MISMATCHES ===
    mismatches = [r for r in results if r.existing_family and r.predicted_family and r.existing_family != r.predicted_family]
    
    if mismatches:
        console.print("\n[bold yellow]═══ MISMATCHES (existing ≠ predicted) ═══[/bold yellow]\n")
        
        for r in mismatches:
            console.print(f"[yellow bold]{r.rcsb_id}.{r.entity_id}[/yellow bold]")
            console.print(f"  Description: {r.description}")
            console.print(f"  Existing: [red]{r.existing_family}[/red] → Predicted: [green]{r.predicted_family}[/green]")
            console.print(f"  Length: {r.sequence_length}, Score: {r.best_score:.1f}, Margin: {r.margin:.1f}")
            scores_str = ", ".join(f"{k}={v:.0f}" for k, v in sorted(r.all_scores.items(), key=lambda x: -x[1]))
            console.print(f"  All scores: {scores_str}")
            console.print()
    
    # === LOW CONFIDENCE TUBULINS ===
    low_conf = [r for r in results if r.confidence in ("weak", "low_margin", "partial") and r.predicted_family]
    
    if low_conf:
        console.print("\n[bold magenta]═══ LOW CONFIDENCE TUBULIN CLASSIFICATIONS ═══[/bold magenta]\n")
        
        for r in sorted(low_conf, key=lambda x: x.margin):
            conf_color = {"weak": "yellow", "low_margin": "yellow", "partial": "magenta"}[r.confidence]
            console.print(f"[{conf_color} bold]{r.rcsb_id}.{r.entity_id}[/{conf_color} bold] [{conf_color}]({r.confidence})[/{conf_color}]")
            console.print(f"  Description: {r.description}")
            console.print(f"  Predicted: {r.predicted_family}, Length: {r.sequence_length}")
            console.print(f"  Score: {r.best_score:.1f}, Margin: {r.margin:.1f}")
            scores_str = ", ".join(f"{k}={v:.0f}" for k, v in sorted(r.all_scores.items(), key=lambda x: -x[1]))
            console.print(f"  All scores: {scores_str}")
            console.print()
    
    # === NO HITS ===
    no_hits = [r for r in results if r.confidence == "no_hit"]
    
    if no_hits:
        console.print("\n[bold red]═══ NO HIT PROTEINS ═══[/bold red]\n")
        console.print(f"These proteins produced no HMM hits at all ({len(no_hits)} total):\n")
        
        table = Table(show_header=True, header_style="bold red")
        table.add_column("PDB", width=6)
        table.add_column("Entity", width=6)
        table.add_column("Length", width=7)
        table.add_column("Existing", width=10)
        table.add_column("Description")
        
        for r in sorted(no_hits, key=lambda x: x.sequence_length):
            existing_str = r.existing_family or "-"
            existing_style = "yellow" if r.existing_family else "dim"
            table.add_row(
                r.rcsb_id,
                r.entity_id,
                str(r.sequence_length),
                f"[{existing_style}]{existing_str}[/{existing_style}]",
                r.description
            )
        
        console.print(table)
    
    # === NOT TUBULIN ===
    not_tubulin = [r for r in results if r.confidence == "not_tubulin"]
    
    console.print("\n[bold cyan]═══ NON-TUBULIN PROTEINS ═══[/bold cyan]\n")
    console.print(f"Proteins with best score < 100 ({len(not_tubulin)} total):\n")
    
    console.print("[bold]By category:[/bold]")
    for cat, count in sorted(stats["non_tubulin_categories"].items(), key=lambda x: -x[1]):
        console.print(f"  {cat:20s}: {count:4d}")
    
    console.print()
    
    # Group by category for display
    by_category = defaultdict(list)
    for r in not_tubulin:
        cat = categorize_non_tubulin(r.description)
        by_category[cat].append(r)
    
    for cat in sorted(by_category.keys(), key=lambda c: -len(by_category[c])):
        cat_results = by_category[cat]
        console.print(f"\n[bold cyan]── {cat} ({len(cat_results)}) ──[/bold cyan]")
        
        table = Table(show_header=True, header_style="bold")
        table.add_column("PDB", width=6)
        table.add_column("Entity", width=6)
        table.add_column("Length", width=7)
        table.add_column("Best Score", width=10)
        table.add_column("Description")
        
        for r in sorted(cat_results, key=lambda x: x.best_score):
            score_style = "green" if r.best_score < 50 else "yellow" if r.best_score < 80 else "red"
            table.add_row(
                r.rcsb_id,
                r.entity_id,
                str(r.sequence_length),
                f"[{score_style}]{r.best_score:.1f}[/{score_style}]",
                r.description[:70]
            )
        
        console.print(table)
    
    # === SCORE DISTRIBUTION ===
    console.print("\n[bold]Non-tubulin score distribution:[/bold]")
    scores = [r.best_score for r in not_tubulin]
    if scores:
        brackets = [(0, 25), (25, 50), (50, 75), (75, 100)]
        for low, high in brackets:
            count = sum(1 for s in scores if low <= s < high)
            bar = "█" * count
            console.print(f"  {low:3d}-{high:3d}: {count:4d} {bar}")
    
    # === HIGH CONFIDENCE TUBULINS (just counts) ===
    high_conf = [r for r in results if r.confidence == "high"]
    console.print(f"\n[bold green]═══ HIGH CONFIDENCE TUBULINS: {len(high_conf)} ═══[/bold green]")
    
    by_family = defaultdict(list)
    for r in high_conf:
        by_family[r.predicted_family].append(r)
    
    for family in ["alpha", "beta", "gamma", "delta", "epsilon"]:
        count = len(by_family.get(family, []))
        console.print(f"  {family:10s}: {count:4d}")
    
    console.print("\n[dim]Use --verbose to see all high-confidence classifications[/dim]")


def main():
    parser = argparse.ArgumentParser(description="PDB Structure HMM Validation")
    subparsers = parser.add_subparsers(dest="command", help="Commands")
    
    # Run command
    run_parser = subparsers.add_parser("run", help="Run analysis and save report")
    run_parser.add_argument("--output", "-o", type=Path, default=DEFAULT_REPORT_PATH,
                           help="Output JSON path")
    
    # Show command
    show_parser = subparsers.add_parser("show", help="Display existing report")
    show_parser.add_argument("--file", "-f", type=Path, default=DEFAULT_REPORT_PATH,
                            help="Report JSON path")
    show_parser.add_argument("--verbose", "-v", action="store_true",
                            help="Show all classifications including high-confidence")
    
    args = parser.parse_args()
    
    if args.command == "run":
        report = run_analysis()
        
        with open(args.output, "w") as f:
            json.dump(report, f, indent=2)
        console.print(f"\n[green]Report saved to {args.output}[/green]\n")
        
        display_report(report)
        
    elif args.command == "show":
        if not args.file.exists():
            console.print(f"[red]Report not found: {args.file}[/red]")
            console.print("Run 'python validate_pdb_structures.py run' first")
            return
        
        with open(args.file) as f:
            report = json.load(f)
        
        display_report(report)
        
    else:
        parser.print_help()

if __name__ == "__main__":
    main()