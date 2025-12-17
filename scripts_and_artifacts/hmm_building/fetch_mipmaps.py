"""
MAP/MIP Sequence Fetcher - Robust Version v2
Uses InterPro domains, excludes fragments, splits large families.
"""
import requests
import time
import os
from dataclasses import dataclass, field


@dataclass
class FamilyConfig:
    desc: str
    interpro: list[str] = field(default_factory=list)
    genes: list[str] = field(default_factory=list)
    min_len: int = 50
    max_len: int = 10000
    reviewed_only: bool = False
    exclude_fragments: bool = True  # NEW: exclude fragments by default
    exclude_terms: list[str] = field(default_factory=list)


MAP_FAMILIES: dict[str, FamilyConfig] = {
    
    # =========================================================================
    # CLASSICAL STRUCTURAL MAPs
    # =========================================================================
    
    # Split tau/MAP2/MAP4 - they have very different sizes
    "tau": FamilyConfig(
        desc="Tau (MAPT)",
        interpro=["IPR027324"],
        genes=["MAPT"],
        min_len=350, max_len=850,  # Tau is ~400-800
    ),
    
    "map2": FamilyConfig(
        desc="MAP2",
        genes=["MAP2"],
        min_len=1500, max_len=2000,  # MAP2 is ~1800
    ),
    
    "map4": FamilyConfig(
        desc="MAP4",
        genes=["MAP4"],
        min_len=1000, max_len=1400,  # MAP4 is ~1100
    ),
    
    # MAP1 - split heavy and light chains
    "map1_heavy": FamilyConfig(
        desc="MAP1A/MAP1B heavy chain",
        genes=["MAP1A", "MAP1B"],
        min_len=2400, max_len=3000,
        exclude_terms=["light chain", "LC3"],
    ),
    
    "map1s": FamilyConfig(
        desc="MAP1S",
        genes=["MAP1S"],
        min_len=1000, max_len=1200,
    ),
    
    "map7": FamilyConfig(
        desc="MAP7/Ensconsin family",
        interpro=["IPR028230"],
        genes=["MAP7", "MAP7D1", "MAP7D2", "MAP7D3"],
        min_len=600, max_len=900,  # tightened
    ),
    
    "doublecortin": FamilyConfig(
        desc="Doublecortin domain family",
        interpro=["IPR003533"],
        genes=["DCX", "DCLK1", "DCLK2", "DCLK3"],
        min_len=350, max_len=750,  # DCX domain proteins, exclude kinase fusions
    ),
    
    # =========================================================================
    # PLUS-END TRACKING PROTEINS (+TIPs)
    # =========================================================================
    
    "eb_family": FamilyConfig(
        desc="EB1/EB3/MAPRE family",
        interpro=["IPR000357", "IPR017975"],
        genes=["MAPRE1", "MAPRE2", "MAPRE3", "EB1", "BIM1", "MAL3"],
        min_len=250, max_len=330,  # EBs are very conserved ~268
    ),
    
    # Split TOG domain proteins
    "ckap5_chtog": FamilyConfig(
        desc="CKAP5/chTOG/XMAP215 (polymerase)",
        genes=["CKAP5"],
        min_len=1800, max_len=2200,  # chTOG is ~2000
    ),
    
    "clasp": FamilyConfig(
        desc="CLASP1/CLASP2",
        genes=["CLASP1", "CLASP2"],
        min_len=1400, max_len=1700,  # CLASPs are ~1500
    ),
    
    # Split CLIP proteins
    "clip170": FamilyConfig(
        desc="CLIP-170/CLIP1",
        genes=["CLIP1"],
        min_len=1400, max_len=1600,
    ),
    
    "clip115": FamilyConfig(
        desc="CLIP-115/CLIP2", 
        genes=["CLIP2"],
        min_len=1000, max_len=1200,
    ),
    
    "tacc": FamilyConfig(
        desc="TACC family",
        interpro=["IPR026939"],
        genes=["TACC1", "TACC2", "TACC3"],
        min_len=700, max_len=900,  # tightened around TACC domain
    ),
    
    # =========================================================================
    # MINUS-END PROTEINS
    # =========================================================================
    
    # Split CAMSAPs by size
    "camsap1": FamilyConfig(
        desc="CAMSAP1",
        genes=["CAMSAP1"],
        min_len=1500, max_len=1700,
    ),
    
    "camsap2": FamilyConfig(
        desc="CAMSAP2",
        genes=["CAMSAP2"],
        min_len=1400, max_len=1600,
    ),
    
    "camsap3": FamilyConfig(
        desc="CAMSAP3/Nezha",
        genes=["CAMSAP3"],
        min_len=1100, max_len=1400,
    ),
    
    # Split gamma-TuRC by size
    "gcp2_3": FamilyConfig(
        desc="GCP2/GCP3 (small TuRC components)",
        genes=["TUBGCP2", "TUBGCP3"],
        min_len=850, max_len=1050,
    ),
    
    "gcp4": FamilyConfig(
        desc="GCP4",
        genes=["TUBGCP4"],
        min_len=600, max_len=800,
    ),
    
    "gcp5_6": FamilyConfig(
        desc="GCP5/GCP6 (large TuRC components)",
        genes=["TUBGCP5", "TUBGCP6"],
        min_len=1000, max_len=1700,
    ),
    
    # =========================================================================
    # SPINDLE/MITOTIC MAPs
    # =========================================================================
    
    "tpx2": FamilyConfig(
        desc="TPX2 spindle assembly factor",
        interpro=["IPR027329"],
        genes=["TPX2"],
        min_len=700, max_len=800,  # tightened
    ),
    
    "numa": FamilyConfig(
        desc="NuMA spindle pole protein",
        interpro=["IPR002743"],
        genes=["NUMA1"],
        min_len=2000, max_len=2300,  # tightened
    ),
    
    # =========================================================================
    # DESTABILIZERS
    # =========================================================================
    
    "stathmin": FamilyConfig(
        desc="Stathmin/Op18 family",
        interpro=["IPR002145"],
        genes=["STMN1", "STMN2", "STMN3", "STMN4"],
        min_len=140, max_len=200,  # very conserved
    ),
    
    "katanin_p60": FamilyConfig(
        desc="Katanin p60 catalytic subunit",
        interpro=["IPR041569"],
        genes=["KATNA1", "KATNAL1", "KATNAL2"],
        min_len=480, max_len=560,
    ),
    
    "spastin": FamilyConfig(
        desc="Spastin MT-severing",
        interpro=["IPR015415"],
        genes=["SPAST"],
        min_len=600, max_len=680,
    ),
    
    "fidgetin": FamilyConfig(
        desc="Fidgetin MT-severing",
        genes=["FIGN", "FIGNL1", "FIGNL2"],
        min_len=700, max_len=780,
        reviewed_only=True,
    ),
    
    "kinesin13": FamilyConfig(
        desc="Kinesin-13 depolymerases",
        interpro=["IPR027129"],
        genes=["KIF2A", "KIF2B", "KIF2C"],
        min_len=700, max_len=800,  # tightened
    ),
    
    # =========================================================================
    # BUNDLERS/CROSSLINKERS
    # =========================================================================
    
    "prc1": FamilyConfig(
        desc="PRC1/Ase1 crosslinker",
        interpro=["IPR008634"],
        genes=["PRC1"],
        min_len=580, max_len=650,
    ),
    
    # =========================================================================
    # TUBULIN CODE - WRITERS
    # =========================================================================
    
    # Split TTLLs by function
    "ttll_glutamylase_short": FamilyConfig(
        desc="TTLL glutamylases (short)",
        genes=["TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7"],
        min_len=300, max_len=500,
    ),
    
    "ttll_glutamylase_long": FamilyConfig(
        desc="TTLL glutamylases (long)",
        genes=["TTLL5", "TTLL6"],
        min_len=1000, max_len=1400,
    ),
    
    "atat1": FamilyConfig(
        desc="Alpha-tubulin acetyltransferase",
        interpro=["IPR026920"],
        genes=["ATAT1"],
        min_len=380, max_len=450,
    ),
    
    # =========================================================================
    # TUBULIN CODE - ERASERS
    # =========================================================================
    
    "ccp_deglutamylase": FamilyConfig(
        desc="CCP/AGBL cytosolic carboxypeptidases",
        interpro=["IPR034859"],
        genes=["AGTPBP1", "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5"],
        min_len=1100, max_len=1350,  # tightened
    ),
    
    "vash_detyrosinase": FamilyConfig(
        desc="Vasohibin tubulin detyrosinases",
        interpro=["IPR029484"],
        genes=["VASH1", "VASH2"],
        min_len=350, max_len=400,  # tightened
    ),
    
    # =========================================================================
    # MIPs / LUMINAL PROTEINS
    # =========================================================================
    
    "fap20_cfap20": FamilyConfig(
        desc="FAP20/CFAP20 inner junction protein",
        interpro=["IPR026621"],
        genes=["CFAP20", "FAP20"],
        min_len=180, max_len=240,  # tightened
    ),
    
    "pacrg": FamilyConfig(
        desc="PACRG inner junction protein",
        interpro=["IPR028138"],
        genes=["PACRG"],
        min_len=300, max_len=380,
    ),
    
    "rib72_efhc": FamilyConfig(
        desc="RIB72/EFHC family (DM10 domain)",
        interpro=["IPR018725"],
        genes=["EFHC1", "EFHC2"],
        min_len=620, max_len=760,  # tightened
    ),
    
    "cfap53": FamilyConfig(
        desc="CFAP53 (FAP53 ortholog)",
        interpro=["IPR029607"],
        genes=["CFAP53", "CCDC11"],
        min_len=380, max_len=480,
    ),
    
    "nme7": FamilyConfig(
        desc="NME7 (FAP67 ortholog - NDPK)",
        genes=["NME7"],
        min_len=370, max_len=420,
    ),
    
    "nme8": FamilyConfig(
        desc="NME8/TXNDC3",
        genes=["NME8", "TXNDC3"],
        min_len=560, max_len=620,
    ),
    
    "spag6": FamilyConfig(
        desc="SPAG6 (PF16 ortholog)",
        interpro=["IPR028152"],
        genes=["SPAG6"],
        min_len=480, max_len=540,
    ),
}


def build_query(config: FamilyConfig) -> str:
    """Build UniProt query from config."""
    parts = []
    
    if config.interpro:
        interpro_parts = [f"xref:interpro-{ipr}" for ipr in config.interpro]
        parts.append("(" + " OR ".join(interpro_parts) + ")")
    
    if config.genes:
        gene_parts = [f"gene:{g}" for g in config.genes]
        parts.append("(" + " OR ".join(gene_parts) + ")")
    
    if not parts:
        raise ValueError("Config must have interpro or genes")
    
    query = " OR ".join(parts)
    
    # Length filter
    query = f"({query}) AND length:[{config.min_len} TO {config.max_len}]"
    
    # Exclude fragments
    if config.exclude_fragments:
        query = f"({query}) AND fragment:false"
    
    # Reviewed only
    if config.reviewed_only:
        query = f"({query}) AND reviewed:true"
    
    return query


@dataclass
class UniProtEntry:
    accession: str
    entry_name: str
    protein_name: str
    gene_name: str | None
    organism: str
    tax_id: int
    sequence: str
    length: int
    reviewed: bool

    def to_fasta(self) -> str:
        prefix = "sp" if self.reviewed else "tr"
        header = (
            f">{prefix}|{self.accession}|{self.entry_name} "
            f"{self.protein_name} OS={self.organism} OX={self.tax_id} "
            f"GN={self.gene_name or 'N/A'}"
        )
        seq_lines = [self.sequence[i:i+60] for i in range(0, len(self.sequence), 60)]
        return header + "\n" + "\n".join(seq_lines)


def fetch_uniprot(query: str, limit: int = 500) -> list[dict]:
    """Fetch from UniProt REST API."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    fields = [
        "accession", "id", "protein_name", "gene_names",
        "organism_name", "organism_id", "sequence", "length", "reviewed"
    ]
    
    params = {
        "query": query,
        "format": "json",
        "fields": ",".join(fields),
        "size": min(limit, 500),
    }
    
    results = []
    
    while len(results) < limit:
        try:
            resp = requests.get(base_url, params=params, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            results.extend(data.get("results", []))
            
            link = resp.headers.get("Link", "")
            if 'rel="next"' in link:
                next_url = link.split(";")[0].strip("<>")
                params = {}
                base_url = next_url
                time.sleep(0.3)
            else:
                break
        except requests.exceptions.RequestException as e:
            print(f"  Error: {e}")
            break
    
    return results[:limit]


def parse_entry(raw: dict) -> UniProtEntry:
    gene_name = None
    if raw.get("genes"):
        first = raw["genes"][0]
        if "geneName" in first:
            gene_name = first["geneName"].get("value")

    protein_name = "Unknown"
    pd = raw.get("proteinDescription", {})
    if "recommendedName" in pd:
        protein_name = pd["recommendedName"]["fullName"]["value"]
    elif pd.get("submissionNames"):
        protein_name = pd["submissionNames"][0]["fullName"]["value"]

    return UniProtEntry(
        accession=raw["primaryAccession"],
        entry_name=raw.get("uniProtkbId", raw["primaryAccession"]),
        protein_name=protein_name,
        gene_name=gene_name,
        organism=raw.get("organism", {}).get("scientificName", "Unknown"),
        tax_id=raw.get("organism", {}).get("taxonId", 0),
        sequence=raw.get("sequence", {}).get("value", ""),
        length=raw.get("sequence", {}).get("length", 0),
        reviewed=raw.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)",
    )


def filter_entries(entries: list[UniProtEntry], config: FamilyConfig) -> list[UniProtEntry]:
    """Filter out entries matching exclude terms."""
    if not config.exclude_terms:
        return entries
    
    filtered = []
    for e in entries:
        text = f"{e.protein_name} {e.gene_name or ''}".lower()
        if not any(term.lower() in text for term in config.exclude_terms):
            filtered.append(e)
    
    return filtered


def fetch_family(key: str, config: FamilyConfig, limit: int = 1000) -> list[UniProtEntry]:
    """Fetch sequences for a MAP/MIP family."""
    query = build_query(config)
    
    print(f"\nFetching {key} ({config.desc})...")
    print(f"  Query: {query[:100]}{'...' if len(query) > 100 else ''}")
    
    raw = fetch_uniprot(query, limit=limit)
    entries = [parse_entry(e) for e in raw]
    
    before = len(entries)
    entries = filter_entries(entries, config)
    if before != len(entries):
        print(f"  Filtered: {before} -> {len(entries)} (exclude_terms)")
    
    reviewed = sum(1 for e in entries if e.reviewed)
    print(f"  Retrieved: {len(entries)} sequences ({reviewed} reviewed)")
    
    if entries:
        lengths = [e.length for e in entries]
        print(f"  Lengths: {min(lengths)}-{max(lengths)} (median {sorted(lengths)[len(lengths)//2]})")
        
        from collections import Counter
        genes = Counter(e.gene_name for e in entries if e.gene_name)
        top_genes = genes.most_common(5)
        print(f"  Top genes: {', '.join(f'{g}:{n}' for g, n in top_genes)}")
    
    return entries


def write_fasta(entries: list[UniProtEntry], path: str):
    with open(path, "w") as f:
        for e in entries:
            f.write(e.to_fasta() + "\n")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Fetch MAP/MIP sequences from UniProt")
    parser.add_argument("--output-dir", default="map_sequences", help="Output directory")
    parser.add_argument("--families", nargs="+", default=None, help="Specific families to fetch")
    parser.add_argument("--limit", type=int, default=1000, help="Max sequences per family")
    parser.add_argument("--list", action="store_true", help="List available families")
    
    args = parser.parse_args()
    
    if args.list:
        print("Available families:\n")
        for key, config in MAP_FAMILIES.items():
            ipr = ", ".join(config.interpro) if config.interpro else "none"
            flags = []
            if config.reviewed_only:
                flags.append("reviewed-only")
            if not config.exclude_fragments:
                flags.append("includes-fragments")
            flags_str = f" [{', '.join(flags)}]" if flags else ""
            print(f"  {key:<25} {config.desc}{flags_str}")
            print(f"    Len: {config.min_len}-{config.max_len}, InterPro: {ipr}")
        exit(0)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    families = args.families or list(MAP_FAMILIES.keys())
    
    print(f"{'='*60}")
    print(f"MAP/MIP Sequence Fetch - {len(families)} families")
    print(f"{'='*60}")
    
    summary = []
    
    for key in families:
        if key not in MAP_FAMILIES:
            print(f"\nWARNING: Unknown family '{key}', skipping")
            continue
        
        config = MAP_FAMILIES[key]
        entries = fetch_family(key, config, limit=args.limit)
        
        if entries:
            out_file = os.path.join(args.output_dir, f"map_{key}.fasta")
            write_fasta(entries, out_file)
            lengths = [e.length for e in entries]
            summary.append((key, len(entries), sum(1 for e in entries if e.reviewed), 
                           min(lengths), max(lengths)))
        else:
            print(f"  WARNING: No results for {key}")
            summary.append((key, 0, 0, 0, 0))
        
        time.sleep(0.5)
    
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"{'Family':<25} {'Total':>6} {'Reviewed':>8} {'Len range':<15}")
    print("-" * 55)
    for key, total, reviewed, min_l, max_l in summary:
        ratio = f"{max_l/min_l:.1f}x" if min_l > 0 else "N/A"
        print(f"{key:<25} {total:>6} {reviewed:>8} {min_l}-{max_l} ({ratio})")
    
    # Flag problematic length ratios
    print("\nLength ratio check (>2x may cause MSA issues):")
    for key, total, reviewed, min_l, max_l in summary:
        if min_l > 0 and max_l / min_l > 2.0:
            print(f"  WARNING: {key} has {max_l/min_l:.1f}x length ratio")
    
    print(f"\nDone. Check '{args.output_dir}/' folder.")