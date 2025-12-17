"""
MAP/MIP Sequence Processing Pipeline
Cluster, balance, and analyze MAP family sequences for HMM building.
"""
import subprocess
import os
import re
import random
from dataclasses import dataclass
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from lib.etl.libtax import Taxid


@dataclass
class ClusteringConfig:
    identity_threshold: float = 0.70
    word_size: int = 5
    min_coverage: float = 0.8
    memory: int = 4000
    threads: int = 4


@dataclass
class BalancingConfig: 

      max_per_phylum  : int       = 5
      max_per_kingdom: int | None = None
      min_total       : int       = 30
      seed            : int       = 42


# MAPs are more divergent than tubulins, so lower thresholds
CLUSTERING_DEFAULTS: dict[str, ClusteringConfig] = {
    "classical_map_tau": ClusteringConfig(identity_threshold=0.70),
    "map1_family"      : ClusteringConfig(identity_threshold=0.70),
    "map7_ensconsin"   : ClusteringConfig(identity_threshold=0.70),
    "map6_stop"        : ClusteringConfig(identity_threshold=0.70),
    "doublecortin"     : ClusteringConfig(identity_threshold=0.70),
    
    # +TIPs - well conserved
    "plus_eb_family"  : ClusteringConfig(identity_threshold=0.75),
    "plus_xmap215_tog": ClusteringConfig(identity_threshold=0.65),   # big, divergent
    "plus_clip"       : ClusteringConfig(identity_threshold=0.70),
    "plus_tacc"       : ClusteringConfig(identity_threshold=0.70),
    
    # Minus-end
    "minus_camsap": ClusteringConfig(identity_threshold=0.65),
    "minus_gamma_turc": ClusteringConfig(identity_threshold=0.70),
    
    # Spindle
    "spindle_tpx2": ClusteringConfig(identity_threshold=0.65),
    "spindle_numa": ClusteringConfig(identity_threshold=0.60),  # very divergent
    "spindle_astrin": ClusteringConfig(identity_threshold=0.65),
    
    # Destabilizers
    "destab_stathmin": ClusteringConfig(identity_threshold=0.75),  # small, conserved
    "destab_severing": ClusteringConfig(identity_threshold=0.65),
    "destab_kinesin13": ClusteringConfig(identity_threshold=0.70),
    
    # Enzymes - catalytic domains conserved
    "code_ttll": ClusteringConfig(identity_threshold=0.65),
    "code_atat": ClusteringConfig(identity_threshold=0.70),
    "code_ccp_deglutamylase": ClusteringConfig(identity_threshold=0.65),
    "code_vash_detyrosinase": ClusteringConfig(identity_threshold=0.70),
    "code_hdac6_deacetylase": ClusteringConfig(identity_threshold=0.70),
    
    # Crosslinkers
    "bundler_prc1": ClusteringConfig(identity_threshold=0.70),
    "spectraplakin": ClusteringConfig(identity_threshold=0.55),  # huge, divergent
}

BALANCING_DEFAULTS = BalancingConfig(max_per_phylum=4, max_per_kingdom=40, min_total=50)


def parse_tax_id_from_header(header: str) -> int | None:
    match = re.search(r'OX=(\d+)', header)
    return int(match.group(1)) if match else None


def parse_gene_from_header(header: str) -> str | None:
    match = re.search(r'GN=(\S+)', header)
    return match.group(1) if match else None


def parse_organism_from_header(header: str) -> str | None:
    match = re.search(r'OS=([^=]+?)\s+OX=', header)
    return match.group(1).strip() if match else None


def get_phylum_for_taxid(tax_id: int) -> str | None:
    if not tax_id:
        return None
    try:
        phylum_id = Taxid.coerce_to_rank(tax_id, "phylum")
        return Taxid.get_name(phylum_id) if phylum_id else None
    except (LookupError, IndexError, ValueError, Exception):
        return None



def get_kingdom_for_taxid(tax_id: int) -> str | None:
    if not tax_id:
        return None
    try:
        return Taxid.superkingdom(tax_id)
    except (LookupError, ValueError, Exception):
        return None


# -----------------------------------------------------------------------------
# Analysis
# -----------------------------------------------------------------------------

def analyze_raw_sequences(fasta_path: str, family: str) -> dict:
    """Analyze a raw FASTA file before any processing."""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    
    if not records:
        print(f"  No sequences found in {fasta_path}")
        return {}
    
    lengths = [len(r.seq) for r in records]
    tax_ids = [parse_tax_id_from_header(r.description) for r in records]
    genes = [parse_gene_from_header(r.description) for r in records]
    organisms = [parse_organism_from_header(r.description) for r in records]
    
    # Kingdom breakdown
    kingdoms = defaultdict(int)
    phyla = defaultdict(int)
    for tid in tax_ids:
        if tid:
            k = get_kingdom_for_taxid(tid)
            kingdoms[k or "unknown"] += 1
            p = get_phylum_for_taxid(tid)
            if p:
                phyla[p] += 1
        else:
            kingdoms["unknown"] += 1
    
    # Gene name breakdown
    gene_counts = defaultdict(int)
    for g in genes:
        gene_counts[g or "N/A"] += 1
    
    # Reviewed vs unreviewed
    reviewed = sum(1 for r in records if r.description.startswith("sp|"))
    
    print(f"\n{'='*60}")
    print(f"RAW ANALYSIS: {family}")
    print(f"{'='*60}")
    print(f"Total sequences: {len(records)}")
    print(f"Reviewed (Swiss-Prot): {reviewed} ({100*reviewed/len(records):.1f}%)")
    
    print(f"\nLength distribution:")
    print(f"  Min: {min(lengths)}, Max: {max(lengths)}, Median: {sorted(lengths)[len(lengths)//2]}")
    print(f"  Mean: {sum(lengths)/len(lengths):.0f}")
    
    print(f"\nKingdom breakdown:")
    for k, n in sorted(kingdoms.items(), key=lambda x: -x[1]):
        print(f"  {k}: {n} ({100*n/len(records):.1f}%)")
    
    print(f"\nTop 10 phyla:")
    for p, n in sorted(phyla.items(), key=lambda x: -x[1])[:10]:
        print(f"  {p}: {n}")
    
    print(f"\nGene names:")
    for g, n in sorted(gene_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {g}: {n}")
    
    print(f"\nTop organisms:")
    org_counts = defaultdict(int)
    for o in organisms:
        org_counts[o or "Unknown"] += 1
    for o, n in sorted(org_counts.items(), key=lambda x: -x[1])[:10]:
        print(f"  {o}: {n}")
    
    return {
        "family": family,
        "total": len(records),
        "reviewed": reviewed,
        "length_min": min(lengths),
        "length_max": max(lengths),
        "length_median": sorted(lengths)[len(lengths)//2],
        "kingdoms": dict(kingdoms),
        "top_phyla": dict(sorted(phyla.items(), key=lambda x: -x[1])[:10]),
        "gene_names": dict(gene_counts),
    }


def analyze_cluster_composition(
    raw_fasta: str,
    clstr_file: str,
    top_n_clusters: int = 5,
):
    """Look inside the largest clusters to understand what's being collapsed."""
    # Load all sequences - handle both sp|ACC|NAME and tr|ACC|NAME formats
    records = {}
    for r in SeqIO.parse(raw_fasta, "fasta"):
        parts = r.id.split("|")
        if len(parts) >= 2:
            acc = parts[1]
        else:
            acc = r.id
        records[acc] = r
    
    # Parse clusters
    clusters = []
    current_cluster = {"id": None, "members": []}
    
    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current_cluster["members"]:
                    clusters.append(current_cluster)
                current_cluster = {"id": line.strip(), "members": []}
            else:
                match = re.search(r'>([^|]+\|)?([A-Z0-9]+)', line)
                if match:
                    acc = match.group(2)
                    is_rep = line.strip().endswith("*")
                    current_cluster["members"].append({"acc": acc, "is_rep": is_rep})
        
        if current_cluster["members"]:
            clusters.append(current_cluster)
    
    clusters.sort(key=lambda c: len(c["members"]), reverse=True)
    
    print(f"\n{'Cluster Analysis':-^60}")
    print(f"Total clusters: {len(clusters)}")
    print(f"Singletons: {sum(1 for c in clusters if len(c['members']) == 1)}")
    
    sizes = [len(c["members"]) for c in clusters]
    print(f"Largest: {max(sizes)}, Mean: {sum(sizes)/len(sizes):.1f}")
    
    print(f"\nTop {top_n_clusters} clusters:\n")
    
    for cluster in clusters[:top_n_clusters]:
        print(f"{cluster['id']} - {len(cluster['members'])} members")
        
        by_organism = defaultdict(list)
        gene_names = set()
        kingdoms = defaultdict(int)
        
        for member in cluster["members"]:
            acc = member["acc"]
            if acc in records:
                rec = records[acc]
                desc = rec.description
                
                organism = parse_organism_from_header(desc) or "Unknown"
                gene = parse_gene_from_header(desc) or "N/A"
                tax_id = parse_tax_id_from_header(desc)
                kingdom = get_kingdom_for_taxid(tax_id) if tax_id else "unknown"
                
                by_organism[organism].append(gene)
                gene_names.add(gene)
                kingdoms[kingdom] += 1
                
                if member["is_rep"]:
                    print(f"  REP: {acc} | {gene} | {organism}")
        
        print(f"  Kingdoms: {dict(kingdoms)}")
        print(f"  Unique organisms: {len(by_organism)}")
        print(f"  Gene names: {gene_names}")
        
        # Top organisms in cluster
        top_orgs = sorted(by_organism.items(), key=lambda x: -len(x[1]))[:3]
        for org, genes in top_orgs:
            print(f"    {org}: {len(genes)} seqs")
        print()


# -----------------------------------------------------------------------------
# Clustering
# -----------------------------------------------------------------------------

def run_cdhit(input_fasta: str, output_fasta: str, config: ClusteringConfig) -> str:
    # Set word size based on identity threshold
    if config.identity_threshold >= 0.7:
        word_size = 5
    elif config.identity_threshold >= 0.6:
        word_size = 4
    elif config.identity_threshold >= 0.5:
        word_size = 3
    else:
        word_size = 2
    
    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(config.identity_threshold),
        "-n", str(word_size),
        "-aL", str(config.min_coverage),
        "-d", "0",
        "-M", str(config.memory),
        "-T", str(config.threads),
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"CD-HIT stderr: {result.stderr}")
        raise RuntimeError(f"CD-HIT failed")
    
    return output_fasta + ".clstr"


# -----------------------------------------------------------------------------
# Balancing
# -----------------------------------------------------------------------------
def balance_or_downsample(
    records: list[SeqRecord],
    max_seqs: int = 300,
    min_seqs: int = 5,
    seed: int = 42,
) -> list[SeqRecord] | None:
    """
    Simple downsampling that prefers reviewed sequences.
    Returns None if fewer than min_seqs available (family should be skipped).
    """
    if len(records) < min_seqs:
        return None
    
    random.seed(seed)
    
    if len(records) <= max_seqs:
        return records
    
    # Prefer reviewed (Swiss-Prot) sequences
    reviewed = [r for r in records if r.description.startswith("sp|")]
    unreviewed = [r for r in records if not r.description.startswith("sp|")]
    
    # Take all reviewed + sample from unreviewed
    if len(reviewed) >= max_seqs:
        return random.sample(reviewed, max_seqs)
    
    n_unreviewed = max_seqs - len(reviewed)
    sampled = reviewed + random.sample(unreviewed, min(n_unreviewed, len(unreviewed)))
    return sampled


def process_map_family(
    input_fasta: str,
    output_dir: str,
    family: str,
    clustering_config: ClusteringConfig | None = None,
    skip_clustering: bool = False,
    analyze_clusters: bool = True,
    min_seqs: int = 5,
    max_seqs: int = 300,
) -> dict | None:
    """Full pipeline for one MAP family. Returns None if family should be skipped."""
    if clustering_config is None:
        clustering_config = CLUSTERING_DEFAULTS.get(family, ClusteringConfig())
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze raw
    stats = analyze_raw_sequences(input_fasta, family)
    if not stats:
        return None
    
    raw_records = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Check minimum before even clustering
    if len(raw_records) < min_seqs:
        print(f"\n  SKIPPING {family}: only {len(raw_records)} raw sequences (need >= {min_seqs})")
        return None
    
    # Cluster
    clustered_path = os.path.join(output_dir, f"map_{family}_nr.fasta")
    
    if skip_clustering and os.path.exists(clustered_path):
        print(f"\nUsing existing clustered file: {clustered_path}")
        clustered_records = list(SeqIO.parse(clustered_path, "fasta"))
        clstr_path = clustered_path + ".clstr"
    else:
        print(f"\nClustering at {clustering_config.identity_threshold*100:.0f}% identity...")
        clstr_path = run_cdhit(input_fasta, clustered_path, clustering_config)
        clustered_records = list(SeqIO.parse(clustered_path, "fasta"))
        
        print(f"  {len(raw_records)} -> {len(clustered_records)} sequences")
        print(f"  Reduction: {100*(1 - len(clustered_records)/len(raw_records)):.1f}%")
    
    # Check minimum after clustering
    if len(clustered_records) < min_seqs:
        print(f"\n  SKIPPING {family}: only {len(clustered_records)} sequences after clustering (need >= {min_seqs})")
        # Clean up the files we just made
        if os.path.exists(clustered_path):
            os.remove(clustered_path)
        if os.path.exists(clstr_path):
            os.remove(clstr_path)
        return None
    
    if analyze_clusters and os.path.exists(clstr_path):
        analyze_cluster_composition(input_fasta, clstr_path, top_n_clusters=3)
    
    # Downsample if needed
    final_records = balance_or_downsample(clustered_records, max_seqs=max_seqs, min_seqs=min_seqs)
    
    if final_records is None:
        print(f"\n  SKIPPING {family}: too few sequences")
        return None
    
    if len(final_records) != len(clustered_records):
        print(f"\nDownsampling: {len(clustered_records)} -> {len(final_records)} sequences")
    
    # Write final
    final_path = os.path.join(output_dir, f"map_{family}_final.fasta")
    with open(final_path, "w") as f:
        SeqIO.write(final_records, f, "fasta")
    print(f"\nWrote {len(final_records)} sequences to {final_path}")
    
    stats["clustered_count"] = len(clustered_records)
    stats["final_count"] = len(final_records)
    
    return stats


def process_all_map_families(
    input_dir: str = "map_sequences",
    output_dir: str = "map_processed",
    families: list[str] | None = None,
    skip_clustering: bool = False,
    min_seqs: int = 5,
    max_seqs: int = 300,
):
    """Process all MAP families."""
    if families is None:
        families = []
        for f in os.listdir(input_dir):
            if f.startswith("map_") and f.endswith(".fasta"):
                family = f.replace("map_", "").replace(".fasta", "")
                families.append(family)
    
    all_stats = []
    skipped = []
    
    for family in sorted(families):
        input_path = os.path.join(input_dir, f"map_{family}.fasta")
        if not os.path.exists(input_path):
            print(f"\nSkipping {family}: file not found")
            continue
        
        stats = process_map_family(
            input_path,
            output_dir,
            family,
            skip_clustering=skip_clustering,
            min_seqs=min_seqs,
            max_seqs=max_seqs,
        )
        if stats:
            all_stats.append(stats)
        else:
            skipped.append(family)
    
    # Summary table
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print(f"{'Family':<25} {'Raw':>6} {'Clust':>6} {'Final':>6} {'Len range':<15}")
    print("-" * 60)
    
    for s in all_stats:
        len_range = f"{s.get('length_min', 0)}-{s.get('length_max', 0)}"
        print(f"{s['family']:<25} {s['total']:>6} {s.get('clustered_count', 0):>6} "
              f"{s.get('final_count', 0):>6} {len_range:<15}")
    
    if skipped:
        print(f"\nSKIPPED ({len(skipped)} families with <{min_seqs} sequences):")
        for fam in skipped:
            print(f"  - {fam}")
    
    return all_stats


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Process MAP/MIP sequences")
    parser.add_argument("--input-dir", default="map_sequences")
    parser.add_argument("--output-dir", default="map_processed")
    parser.add_argument("--families", nargs="+", default=None)
    parser.add_argument("--skip-clustering", action="store_true")
    parser.add_argument("--min-seqs", type=int, default=5, help="Minimum sequences to keep family")
    parser.add_argument("--max-seqs", type=int, default=300, help="Max sequences after downsampling")
    parser.add_argument("--analyze-only", action="store_true")
    
    args = parser.parse_args()
    
    if args.analyze_only:
        for f in os.listdir(args.input_dir):
            if f.startswith("map_") and f.endswith(".fasta"):
                family = f.replace("map_", "").replace(".fasta", "")
                analyze_raw_sequences(os.path.join(args.input_dir, f), family)
    else:
        process_all_map_families(
            args.input_dir,
            args.output_dir,
            args.families,
            args.skip_clustering,
            min_seqs=args.min_seqs,
            max_seqs=args.max_seqs,
        )