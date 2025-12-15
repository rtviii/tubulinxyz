import subprocess
import os
import re
import random
from dataclasses import dataclass
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from lib.etl.libtax import Taxid, get_ncbi


@dataclass
class ClusteringConfig:
    """Configuration for sequence clustering."""
    identity_threshold: float = 0.90
    word_size: int = 5
    min_coverage: float = 0.8
    description_length: int = 0
    memory: int = 4000
    threads: int = 4


@dataclass
class BalancingConfig:
    """Configuration for taxonomic balancing."""
    max_per_phylum: int = 5
    max_per_kingdom: int | None = None  # None = no cap
    min_total: int = 30
    seed: int = 42


# Defaults per family
CLUSTERING_DEFAULTS: dict[str, ClusteringConfig] = {
    "alpha": ClusteringConfig(identity_threshold=0.90),
    "beta": ClusteringConfig(identity_threshold=0.90),
    "gamma": ClusteringConfig(identity_threshold=0.85),
    "delta": ClusteringConfig(identity_threshold=0.80),
    "epsilon": ClusteringConfig(identity_threshold=0.80),
    "zeta": ClusteringConfig(identity_threshold=0.75),
}

BALANCING_DEFAULTS: dict[str, BalancingConfig] = {
    "alpha": BalancingConfig(max_per_phylum=5, max_per_kingdom=50, min_total=80),
    "beta": BalancingConfig(max_per_phylum=5, max_per_kingdom=50, min_total=80),
    "gamma": BalancingConfig(max_per_phylum=5, max_per_kingdom=40, min_total=60),
    "delta": BalancingConfig(max_per_phylum=6, max_per_kingdom=60, min_total=80),
    "epsilon": BalancingConfig(max_per_phylum=3, max_per_kingdom=15, min_total=20),
    "zeta": BalancingConfig(max_per_phylum=3, max_per_kingdom=None, min_total=15),
}


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def parse_tax_id_from_header(header: str) -> int | None:
    """Extract taxonomic ID from FASTA header (OX=XXXXX pattern)."""
    match = re.search(r'OX=(\d+)', header)
    return int(match.group(1)) if match else None


def get_phylum_for_record(record: SeqRecord) -> str | None:
    """Get phylum name from a sequence record."""
    tax_id = parse_tax_id_from_header(record.description)
    if not tax_id:
        return None
    try:
        phylum_id = Taxid.coerce_to_rank(tax_id, "phylum")
        return Taxid.get_name(phylum_id) if phylum_id else None
    except (LookupError, IndexError):
        return None


def get_kingdom_for_record(record: SeqRecord) -> str | None:
    """Get kingdom (superkingdom) from a sequence record."""
    tax_id = parse_tax_id_from_header(record.description)
    if not tax_id:
        return None
    try:
        return Taxid.superkingdom(tax_id)
    except (LookupError, ValueError):
        return None


def get_taxonomic_breakdown(tax_ids: list[int]) -> dict[str, int]:
    """Get kingdom-level breakdown of taxonomic IDs."""
    breakdown = defaultdict(int)
    for tax_id in tax_ids:
        if tax_id is None:
            breakdown["unknown"] += 1
            continue
        try:
            kingdom = Taxid.superkingdom(tax_id)
            breakdown[kingdom] += 1
        except (LookupError, ValueError):
            breakdown["unknown"] += 1
    return dict(breakdown)


def get_phylum_breakdown(tax_ids: list[int], top_n: int = 10) -> dict[str, int]:
    """Get phylum-level breakdown, returns top N phyla."""
    breakdown = defaultdict(int)
    for tax_id in tax_ids:
        if tax_id is None:
            continue
        try:
            phylum_id = Taxid.coerce_to_rank(tax_id, "phylum")
            if phylum_id:
                phylum_name = Taxid.get_name(phylum_id)
                breakdown[phylum_name] += 1
        except (LookupError, IndexError, ValueError):
            pass
    sorted_phyla = sorted(breakdown.items(), key=lambda x: -x[1])
    return dict(sorted_phyla[:top_n])


def load_fasta_with_taxids(fasta_path: str) -> tuple[list[SeqRecord], list[int]]:
    """Load FASTA file and extract tax IDs from headers."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    tax_ids = [parse_tax_id_from_header(r.description) for r in records]
    return records, tax_ids


# -----------------------------------------------------------------------------
# Clustering (CD-HIT)
# -----------------------------------------------------------------------------

def run_cdhit(input_fasta: str, output_fasta: str, config: ClusteringConfig) -> str:
    """Run CD-HIT clustering. Returns path to cluster file."""
    if config.word_size is None:
        if config.identity_threshold >= 0.7:
            word_size = 5
        elif config.identity_threshold >= 0.6:
            word_size = 4
        elif config.identity_threshold >= 0.5:
            word_size = 3
        else:
            word_size = 2
    else:
        word_size = config.word_size

    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(config.identity_threshold),
        "-n", str(word_size),
        "-aL", str(config.min_coverage),
        "-d", str(config.description_length),
        "-M", str(config.memory),
        "-T", str(config.threads),
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"CD-HIT stderr: {result.stderr}")
        raise RuntimeError(f"CD-HIT failed with code {result.returncode}")

    return output_fasta + ".clstr"


def parse_cdhit_clusters(clstr_path: str) -> dict[str, list[str]]:
    """Parse CD-HIT .clstr file."""
    clusters = {}
    current_cluster = []
    current_rep = None

    with open(clstr_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                if current_rep and current_cluster:
                    clusters[current_rep] = current_cluster
                current_cluster = []
                current_rep = None
            else:
                match = re.search(r'>([^\s.]+)', line)
                if match:
                    member_id = match.group(1)
                    current_cluster.append(member_id)
                    if line.endswith("*"):
                        current_rep = member_id

        if current_rep and current_cluster:
            clusters[current_rep] = current_cluster

    return clusters


# -----------------------------------------------------------------------------
# Taxonomic Balancing
# -----------------------------------------------------------------------------

def balance_by_taxonomy(
    records: list[SeqRecord],
    config: BalancingConfig,
) -> list[SeqRecord]:
    """
    Balance sequences by taxonomic representation.
    
    Strategy:
    1. Group sequences by phylum
    2. Take up to max_per_phylum from each, preferring diverse sequences
    3. Cap per kingdom if specified
    4. Ensure minimum total sequences
    """
    random.seed(config.seed)

    # Group by taxonomy
    by_kingdom: dict[str, list[SeqRecord]] = defaultdict(list)
    by_phylum: dict[tuple[str, str], list[SeqRecord]] = defaultdict(list)
    unknown: list[SeqRecord] = []

    for record in records:
        kingdom = get_kingdom_for_record(record)
        phylum = get_phylum_for_record(record)

        if kingdom is None:
            unknown.append(record)
            continue

        by_kingdom[kingdom].append(record)
        phylum_key = phylum or f"{kingdom}_unknown"
        by_phylum[(kingdom, phylum_key)].append(record)

    # Select from each phylum, smallest phyla first (ensures rare groups get represented)
    selected: list[SeqRecord] = []
    kingdom_counts: dict[str, int] = defaultdict(int)
    selected_ids: set[str] = set()

    sorted_phyla = sorted(by_phylum.items(), key=lambda x: len(x[1]))

    for (kingdom, phylum), phylum_records in sorted_phyla:
        # Check kingdom cap
        if config.max_per_kingdom and kingdom_counts[kingdom] >= config.max_per_kingdom:
            continue

        # How many can we take from this phylum?
        available = config.max_per_phylum
        if config.max_per_kingdom:
            available = min(available, config.max_per_kingdom - kingdom_counts[kingdom])

        # Select (randomly for now, could be diversity-based)
        n_take = min(available, len(phylum_records))
        sampled = random.sample(phylum_records, n_take)

        for rec in sampled:
            selected.append(rec)
            selected_ids.add(rec.id)
            kingdom_counts[kingdom] += 1

    # Add unknowns
    selected.extend(unknown)

    # If below minimum, add more from underrepresented groups
    if len(selected) < config.min_total:
        remaining_needed = config.min_total - len(selected)

        # Prefer non-Chordata eukaryotes, then anything else
        for kingdom in by_kingdom:
            if remaining_needed <= 0:
                break

            candidates = [r for r in by_kingdom[kingdom] if r.id not in selected_ids]

            # Deprioritize Chordata
            if kingdom == "eukaryota":
                non_chordata = [r for r in candidates if get_phylum_for_record(r) != "Chordata"]
                if non_chordata:
                    candidates = non_chordata

            n_add = min(remaining_needed, len(candidates))
            if n_add > 0:
                extra = random.sample(candidates, n_add)
                selected.extend(extra)
                selected_ids.update(r.id for r in extra)
                remaining_needed -= n_add

    return selected


# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------

def print_stage_comparison(
    before_records: list[SeqRecord],
    after_records: list[SeqRecord],
    stage_name: str,
    family: str,
):
    """Print before/after comparison for a processing stage."""
    before_tax_ids = [parse_tax_id_from_header(r.description) for r in before_records]
    after_tax_ids = [parse_tax_id_from_header(r.description) for r in after_records]

    before_kingdoms = get_taxonomic_breakdown(before_tax_ids)
    after_kingdoms = get_taxonomic_breakdown(after_tax_ids)

    before_phyla = get_phylum_breakdown(before_tax_ids, top_n=8)
    after_phyla = get_phylum_breakdown(after_tax_ids, top_n=8)

    print(f"\n{'='*70}")
    print(f"{stage_name}: {family.upper()}-TUBULIN")
    print(f"{'='*70}")

    print(f"\nSequence counts:")
    print(f"  Before: {len(before_records)}")
    print(f"  After:  {len(after_records)}")
    reduction = 100 * (1 - len(after_records) / len(before_records)) if before_records else 0
    print(f"  Reduction: {len(before_records) - len(after_records)} ({reduction:.1f}%)")

    print(f"\n{'Kingdom distribution':-^50}")
    print(f"  {'Kingdom':<15} {'Before':>10} {'After':>10} {'Retained':>10}")
    print(f"  {'-'*45}")

    all_kingdoms = set(before_kingdoms.keys()) | set(after_kingdoms.keys())
    for kingdom in sorted(all_kingdoms):
        before_n = before_kingdoms.get(kingdom, 0)
        after_n = after_kingdoms.get(kingdom, 0)
        retained = f"{100*after_n/before_n:.0f}%" if before_n > 0 else "N/A"
        print(f"  {kingdom:<15} {before_n:>10} {after_n:>10} {retained:>10}")

    print(f"\n{'Top phyla (after)':-^50}")
    for phylum, count in after_phyla.items():
        print(f"  {phylum:<35} {count:>5}")


def print_cluster_stats(clstr_path: str):
    """Print cluster statistics."""
    clusters = parse_cdhit_clusters(clstr_path)
    cluster_sizes = [len(members) for members in clusters.values()]

    print(f"\n{'Cluster statistics':-^50}")
    print(f"  Total clusters: {len(clusters)}")
    print(f"  Singletons: {sum(1 for s in cluster_sizes if s == 1)}")
    print(f"  Largest cluster: {max(cluster_sizes)} sequences")
    print(f"  Average cluster size: {sum(cluster_sizes)/len(cluster_sizes):.1f}")


# -----------------------------------------------------------------------------
# Main Pipeline
# -----------------------------------------------------------------------------

def process_tubulin_family(
    input_fasta: str,
    output_dir: str,
    family: str,
    clustering_config: ClusteringConfig | None = None,
    balancing_config: BalancingConfig | None = None,
    skip_clustering: bool = False,
    skip_balancing: bool = False,
) -> dict:
    """
    Full pipeline: cluster -> balance -> output.
    
    Outputs:
        - tubulin_{family}_nr.fasta (clustered only)
        - tubulin_{family}_final.fasta (clustered + balanced)
    """
    if clustering_config is None:
        clustering_config = CLUSTERING_DEFAULTS.get(family, ClusteringConfig())
    if balancing_config is None:
        balancing_config = BALANCING_DEFAULTS.get(family, BalancingConfig())

    os.makedirs(output_dir, exist_ok=True)

    # Load raw
    raw_records, _ = load_fasta_with_taxids(input_fasta)
    print(f"\nLoaded {len(raw_records)} raw sequences for {family}-tubulin")

    # Step 1: Clustering
    clustered_path = os.path.join(output_dir, f"tubulin_{family}_nr.fasta")

    if skip_clustering and os.path.exists(clustered_path):
        print(f"Skipping clustering, loading from {clustered_path}")
        clustered_records, _ = load_fasta_with_taxids(clustered_path)
    else:
        print(f"\nClustering at {clustering_config.identity_threshold*100:.0f}% identity...")
        clstr_path = run_cdhit(input_fasta, clustered_path, clustering_config)
        clustered_records, _ = load_fasta_with_taxids(clustered_path)
        print_stage_comparison(raw_records, clustered_records, "CLUSTERING", family)
        print_cluster_stats(clstr_path)

    # Step 2: Taxonomic balancing
    final_path = os.path.join(output_dir, f"tubulin_{family}_final.fasta")

    if skip_balancing:
        print("Skipping balancing")
        final_records = clustered_records
    else:
        print(f"\nBalancing (max {balancing_config.max_per_phylum}/phylum, "
              f"max {balancing_config.max_per_kingdom}/kingdom)...")
        final_records = balance_by_taxonomy(clustered_records, balancing_config)
        print_stage_comparison(clustered_records, final_records, "BALANCING", family)

    # Write final output
    with open(final_path, "w") as f:
        SeqIO.write(final_records, f, "fasta")
    print(f"\nWrote {len(final_records)} sequences to {final_path}")

    # Summary stats
    final_tax_ids = [parse_tax_id_from_header(r.description) for r in final_records]

    return {
        "family": family,
        "raw_count": len(raw_records),
        "clustered_count": len(clustered_records),
        "final_count": len(final_records),
        "clustering_threshold": clustering_config.identity_threshold,
        "kingdoms": get_taxonomic_breakdown(final_tax_ids),
        "top_phyla": get_phylum_breakdown(final_tax_ids, top_n=5),
    }


def process_all_families(
    input_dir: str,
    output_dir: str,
    families: list[str] | None = None,
    skip_clustering: bool = False,
    skip_balancing: bool = False,
):
    """Process all tubulin families."""
    if families is None:
        families = ["alpha", "beta", "gamma", "delta", "epsilon", "zeta"]

    all_stats = []

    for family in families:
        input_path = os.path.join(input_dir, f"tubulin_{family}_raw.fasta")

        if not os.path.exists(input_path):
            print(f"\nSkipping {family}: {input_path} not found")
            continue

        stats = process_tubulin_family(
            input_path,
            output_dir,
            family,
            skip_clustering=skip_clustering,
            skip_balancing=skip_balancing,
        )
        all_stats.append(stats)

    # Final summary
    print(f"\n{'='*70}")
    print("FINAL SUMMARY")
    print(f"{'='*70}")
    print(f"{'Family':<10} {'Raw':>8} {'Clustered':>10} {'Final':>8} {'Kingdoms breakdown':<30}")
    print("-" * 70)

    for s in all_stats:
        kingdoms_str = ", ".join(f"{k}:{v}" for k, v in sorted(s['kingdoms'].items()))
        print(f"{s['family']:<10} {s['raw_count']:>8} {s['clustered_count']:>10} "
              f"{s['final_count']:>8} {kingdoms_str:<30}")

    return all_stats


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Cluster and balance tubulin FASTA files")
    parser.add_argument("--input-dir", default="./raw_fastas", help="Directory with raw FASTA files")
    parser.add_argument("--output-dir", default="./processed", help="Output directory")
    parser.add_argument("--families", nargs="+", default=None, help="Families to process")
    parser.add_argument("--skip-clustering", action="store_true", help="Skip clustering step")
    parser.add_argument("--skip-balancing", action="store_true", help="Skip balancing step")

    args = parser.parse_args()

    process_all_families(
        args.input_dir,
        args.output_dir,
        args.families,
        skip_clustering=args.skip_clustering,
        skip_balancing=args.skip_balancing,
    )