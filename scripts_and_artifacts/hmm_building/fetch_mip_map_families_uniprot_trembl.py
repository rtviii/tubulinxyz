"""
MAP/MIP Sequence Fetcher
Adapts the user's Tubulin fetcher to handle diverse MAP families.
Uses Hybrid Queries (Gene Names OR Pfam Domains) to maximize coverage.
"""
import requests
import time
import json
from dataclasses import dataclass
from typing import Iterator, Dict, Any
import re
import os

# --- Configuration: The "Master List" of MAP Families ---

# --- Configuration: The "Master List" of MAP Families ---

def build_query(gene_names: list[str], interpro_id: str = None, extra_terms: str = None) -> str:
    """
    Constructs a UniProt query: (Genes) OR (InterPro)
    FIX: Uses 'interpro:IPRxxxxx' which is the correct field search syntax.
    """
    parts = []
    
    # Gene name part: (gene:A OR gene:B ...)
    if gene_names:
        genes_str = " OR ".join([f"gene:{g}" for g in gene_names])
        parts.append(f"({genes_str})")
        
    # InterPro ID part: interpro:IPRxxxxx
    if interpro_id:
        # The correct field name in UniProt API is 'interpro'
        parts.append(f"(interpro:{interpro_id})")
        
    # Combine with OR to catch either specific names OR structural homologs
    base_query = " OR ".join(parts)
    
    if extra_terms:
        return f"({base_query}) AND {extra_terms}"
    return base_query

# The Definitions (Updated to use plain InterPro IDs)
MAP_FAMILIES = {
    # --- Level 1: Lattice Binders ---
    "classical_map_tau": {
        "query": build_query(["MAPT", "MAP2", "MAP4"], "IPR001084"), 
        "min_len": 200, "max_len": 3000, 
        "desc": "Tau/MAP2/MAP4 family"
    },
    "doublecortin": {
        "query": build_query(["DCX", "DCLK1", "DCLK2"], "IPR003533"), 
        "min_len": 300, "max_len": 1500,
        "desc": "Doublecortin family"
    },
    "bundler_prc1": {
        "query": build_query(["PRC1", "ASE1", "MAP65"], "IPR000078"), 
        "min_len": 400, "max_len": 1200,
        "desc": "PRC1/MAP65 crosslinkers"
    },

    # --- Level 2: Plus-End Tracking (+TIPs) ---
    "plus_eb_family": {
        "query": build_query(["MAPRE1", "MAPRE2", "MAPRE3", "EB1", "BIM1"], "IPR000628"), 
        "min_len": 150, "max_len": 600,
        "desc": "EB1/EB3 family (The core +TIPs)"
    },
    "plus_xmap215_tog": {
        "query": build_query(["CKAP5", "XMAP215", "CLASP1", "CLASP2"], "IPR002476"), 
        "min_len": 800, "max_len": 5000, 
        "desc": "Polymerases (TOG domain proteins)"
    },
    "plus_clip": {
        "query": build_query(["CLIP1", "CLIP2", "CLIP3", "CLIP4"], "IPR001332"), 
        "min_len": 300, "max_len": 2000,
        "desc": "CLIP-170 family"
    },

    # --- Level 3: Minus-End ---
    "minus_camsap": {
        "query": build_query(["CAMSAP1", "CAMSAP2", "CAMSAP3", "KANSL1"], "IPR029471"), 
        "min_len": 400, "max_len": 2500,
        "desc": "Patronin/CAMSAP family"
    },
    "minus_gamma_turc": {
        "query": build_query(["TUBGCP2", "TUBGCP3", "TUBGCP4", "TUBGCP5", "TUBGCP6"], "IPR010688"), 
        "min_len": 400, "max_len": 2000,
        "desc": "Gamma-tubulin complex proteins (GCPs)"
    },

    # --- Level 4: Destabilizers ---
    "destab_stathmin": {
        "query": build_query(["STMN1", "STMN2", "SCG10"], "IPR000743"), 
        "min_len": 100, "max_len": 400, 
        "desc": "Sequestering Stathmin family"
    },
    "destab_severing_aaa": {
        "query": build_query(["KATNA1", "SPAST", "FIGN"], "IPR003959"), 
        "min_len": 300, "max_len": 1500,
        "desc": "Severing enzymes (Katanin/Spastin)"
    },
    "destab_kinesin13": {
        # Kinesin-13 specific
        "query": '(family:"kinesin-13") AND (interpro:IPR001752)', 
        "min_len": 400, "max_len": 1200,
        "desc": "MCAK/Kinesin-13 depolymerases"
    },

    # --- Level 5: The "Code" (Enzymes) ---
    "code_glutamylase_ttll": {
        "query": build_query(["TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6"], "IPR023253"), 
        "min_len": 300, "max_len": 1500,
        "desc": "Polyglutamylases (TTLL)"
    },
    "code_acetylase_atat": {
        "query": build_query(["ATAT1", "MEC17"], "IPR000182"), 
        "min_len": 150, "max_len": 600,
        "desc": "Alpha-tubulin acetyltransferase"
    }
}




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
    annotation_score: int      # 1-5, 5 is best
    protein_existence: int     # 1-5, 1 is best
    reviewed: bool
    lineage: list[str]         
    
    def to_fasta_header(self) -> str:
        prefix = "sp" if self.reviewed else "tr"
        return (
            f">{prefix}|{self.accession}|{self.entry_name} "
            f"{self.protein_name} "
            f"OS={self.organism} "
            f"OX={self.tax_id} "
            f"GN={self.gene_name or 'N/A'} "
            f"PE={self.protein_existence} "
            f"AS={self.annotation_score} "
            f"LEN={self.length}"
        )
    
    def to_fasta(self) -> str:
        header = self.to_fasta_header()
        seq_lines = [self.sequence[i:i+60] for i in range(0, len(self.sequence), 60)]
        return header + "\n" + "\n".join(seq_lines)

class QueryLogger:
    """Tracks exactly what we queried for reproducibility."""
    def __init__(self, filepath="query_manifest.txt"):
        self.filepath = filepath
        # Initialize file with header
        with open(self.filepath, "w") as f:
            f.write(f"Run Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("Filename | Family_Key | Exact_Query_String | Filters\n")
            f.write("-" * 80 + "\n")
    
    def log(self, filename, family_key, query_str, filters):
        with open(self.filepath, "a") as f:
            f.write(f"{filename} | {family_key} | {query_str} | {filters}\n")


def fetch_uniprot_json(query: str, limit: int = 500) -> list[dict]:
    """Fetch from UniProt REST API with full metadata."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    fields = [
        "accession", "id", "protein_name", "gene_names",
        "organism_name", "organism_id", "sequence", "length",
        "annotation_score", "protein_existence", "reviewed", "lineage"
    ]
    
    params = {
        "query": query,
        "format": "json",
        "fields": ",".join(fields),
        "size": min(limit, 500),
    }
    
    results = []
    
    while True:
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            results.extend(data.get("results", []))
            
            if len(results) >= limit:
                break

            # Handle pagination
            link_header = response.headers.get("Link", "")
            next_match = re.search(r'<([^>]+)>; rel="next"', link_header)
            
            if next_match:
                base_url = next_match.group(1)
                params = {} 
                time.sleep(0.5) 
            else:
                break
        except requests.exceptions.RequestException as e:
            print(f"Error fetching data: {e}")
            break
            
    return results[:limit]


def parse_uniprot_entry(raw: dict) -> UniProtEntry:
    """Parse a UniProt JSON entry into our dataclass."""
    
    gene_name = None
    if "genes" in raw and raw["genes"]:
        first_gene = raw["genes"][0]
        if "geneName" in first_gene:
            gene_name = first_gene["geneName"].get("value")
    
    protein_name = "Unknown"
    if "proteinDescription" in raw:
        pd = raw["proteinDescription"]
        if "recommendedName" in pd:
            protein_name = pd["recommendedName"]["fullName"]["value"]
        elif "submissionNames" in pd and pd["submissionNames"]:
            protein_name = pd["submissionNames"][0]["fullName"]["value"]
            
    lineage = []
    if "lineages" in raw.get("organism", {}):
        lineage = raw["organism"]["lineages"]
    
    return UniProtEntry(
        accession=raw["primaryAccession"],
        entry_name=raw.get("uniProtkbId", raw["primaryAccession"]),
        protein_name=protein_name,
        gene_name=gene_name,
        organism=raw.get("organism", {}).get("scientificName", "Unknown"),
        tax_id=raw.get("organism", {}).get("taxonId", 0),
        sequence=raw.get("sequence", {}).get("value", ""),
        length=raw.get("sequence", {}).get("length", 0),
        annotation_score=raw.get("annotationScore", 0),
        protein_existence=int(raw.get("proteinExistence", "5_uncertain")[-1:]),
        reviewed=raw.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)",
        lineage=lineage,
    )

def fetch_map_family(
    key: str,
    config: dict,
    logger: QueryLogger,
    limit: int = 500
) -> list[UniProtEntry]:
    
    # 1. Build Base Query
    base_query = config["query"]
    
    # 2. Add Length Filters (Critical for MAPs to remove fragments/fusion artifacts)
    len_filter = f" AND length:[{config['min_len']} TO {config['max_len']}]"
    
    # 3. Add Annotation Score Filter (Keep it somewhat high to avoid noise)
    score_filter = " AND annotation_score:[3 TO 5]"
    
    final_query = f"({base_query}){len_filter}{score_filter}"
    
    print(f"Fetching {key} ({config['desc']})...")
    # print(f"  Query: {final_query}") # Uncomment for verbose debug
    
    raw_entries = fetch_uniprot_json(final_query, limit=limit)
    entries = [parse_uniprot_entry(e) for e in raw_entries]
    
    # Log what we did
    filename = f"map_{key}.fasta"
    logger.log(filename, key, base_query, f"len:{config['min_len']}-{config['max_len']}, score>=3")
    
    return entries

def write_fasta(entries: list[UniProtEntry], outpath: str):
    with open(outpath, "w") as f:
        for entry in entries:
            f.write(entry.to_fasta() + "\n")

# --- Main Execution ---

if __name__ == "__main__":
    
    # Initialize Logger
    logger = QueryLogger("map_query_manifest.txt")
    
    # Create output directory
    os.makedirs("map_sequences", exist_ok=True)
    
    print(f"{'='*60}")
    print(f"Starting MAP/MIP Sequence Fetch")
    print(f"Targeting {len(MAP_FAMILIES)} distinct protein families.")
    print(f"{'='*60}\n")
    
    for key, config in MAP_FAMILIES.items():
        
        entries = fetch_map_family(key, config, logger, limit=200) # Limit 200 is good for initial MSAs
        
        if entries:
            out_file = os.path.join("map_sequences", f"map_{key}.fasta")
            write_fasta(entries, out_file)
            
            # Brief Stats
            reviewed_count = sum(1 for e in entries if e.reviewed)
            print(f"  -> Saved {len(entries)} sequences to {out_file}")
            print(f"     ({reviewed_count} Reviewed, {len(entries)-reviewed_count} Unreviewed)")
        else:
            print(f"  -> WARNING: No results found for {key}")
            
        print("-" * 40)
        
        # Respect API rate limits
        time.sleep(1)

    print("\nDone. Check 'map_sequences/' folder and 'map_query_manifest.txt' for details.")