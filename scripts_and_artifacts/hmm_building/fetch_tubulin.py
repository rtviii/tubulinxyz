"""
What this does:
>sp|Q9ZSW1|TBB1_CYAPA Tubulin beta-1 chain OS=Cyanophora paradoxa OX=2762 GN=TUBB1 PE=2 AS=4 LEN=452
"""
import requests
import time
from dataclasses import dataclass
from typing import Iterator
import re

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
    protein_existence: int     # 1-5, 1 is best (confusingly opposite)
    reviewed: bool
    lineage: list[str]         # taxonomic lineage
    
    def to_fasta_header(self) -> str:
        """Custom header with the info we care about."""
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
        # Wrap sequence at 60 chars
        seq_lines = [self.sequence[i:i+60] for i in range(0, len(self.sequence), 60)]
        return header + "\n" + "\n".join(seq_lines)


def fetch_uniprot_json(query: str, limit: int = 500) -> list[dict]:
    """
    Fetch from UniProt REST API with full metadata.
    Returns list of entry dicts.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Fields we want - see https://www.uniprot.org/help/return_fields
    fields = [
        "accession",
        "id",                    # entry name
        "protein_name",
        "gene_names",
        "organism_name",
        "organism_id",
        "sequence",
        "length",
        "annotation_score",
        "protein_existence",
        "reviewed",
        "lineage",               # full taxonomic lineage
    ]
    
    params = {
        "query": query,
        "format": "json",
        "fields": ",".join(fields),
        "size": min(limit, 500),  # API max is 500 per request
    }
    
    results = []
    
    while True:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        
        results.extend(data.get("results", []))
        
        # Handle pagination
        link_header = response.headers.get("Link", "")
        next_match = re.search(r'<([^>]+)>; rel="next"', link_header)
        
        if next_match and len(results) < limit:
            # Follow pagination link
            base_url = next_match.group(1)
            params = {}  # URL already contains params
            time.sleep(0.5)  # Be nice to the API
        else:
            break
    
    return results[:limit]


def parse_uniprot_entry(raw: dict) -> UniProtEntry:
    """Parse a UniProt JSON entry into our dataclass."""
    
    # Gene name extraction (can be complex structure)
    gene_name = None
    if "genes" in raw and raw["genes"]:
        first_gene = raw["genes"][0]
        if "geneName" in first_gene:
            gene_name = first_gene["geneName"].get("value")
    
    # Protein name extraction
    protein_name = "Unknown"
    if "proteinDescription" in raw:
        pd = raw["proteinDescription"]
        if "recommendedName" in pd:
            protein_name = pd["recommendedName"]["fullName"]["value"]
        elif "submissionNames" in pd and pd["submissionNames"]:
            protein_name = pd["submissionNames"][0]["fullName"]["value"]
    
    # Lineage
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
        protein_existence=raw.get("proteinExistence", "5_uncertain")[-1:],  # Extract number
        reviewed=raw.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)",
        lineage=lineage,
    )


# Tubulin-specific queries
# Using InterPro families for more comprehensive results

TUBULIN_QUERIES = {
    # InterPro families - these are well-curated
    "alpha"  : "(xref:interpro-IPR000217) OR (family:tubulin AND protein_name:alpha)",
    "beta"   : "(xref:interpro-IPR000217) OR (family:tubulin AND protein_name:beta)",
    "gamma"  : "(xref:interpro-IPR002452) OR (family:tubulin AND protein_name:gamma)",
    "delta"  : "(family:tubulin AND protein_name:delta)",
    "epsilon": "(family:tubulin AND protein_name:epsilon)",
    "zeta"   : "(family:tubulin AND protein_name:zeta)",
}

# Note: IPR000217 is "Tubulin" (covers alpha/beta), IPR002452 is "Tubulin/FtsZ, GTPase domain"
# The InterPro queries will catch more but may need filtering

# More specific alternative using protein family
TUBULIN_QUERIES_STRICT = {
    "alpha"  : '(family:"tubulin family") AND (protein_name:"alpha" OR gene:TUBA*)',
    "beta"   : '(family:"tubulin family") AND (protein_name:"beta" OR gene:TUBB*)',
    "gamma"  : '(family:"tubulin family") AND (protein_name:"gamma" OR gene:TUBG*)',
    "delta"  : '(family:"tubulin family") AND (protein_name:"delta" OR gene:TUBD*)',
    "epsilon": '(family:"tubulin family") AND (protein_name:"epsilon" OR gene:TUBE*)',
    "zeta"   : '(family:"tubulin family") AND (protein_name:"zeta" OR gene:TUBZ*)',
}


def fetch_tubulin_family(
  family               : str,
  reviewed_only        : bool = False,   # Changed default - TrEMBL has more coverage
  min_annotation_score : int = 3,
  min_length           : int = 400,      # Full tubulin is ~450 aa
  max_length           : int = 600,
  limit                : int = 1000,
) -> list[UniProtEntry]: 
    """
    Fetch tubulin sequences with quality filters.
    """
    query = TUBULIN_QUERIES_STRICT.get(family)
    if not query:
        raise ValueError(f"Unknown family: {family}")
    
    # Add filters to query
    if reviewed_only:
        query += " AND reviewed:true"
    
    query += f" AND length:[{min_length} TO {max_length}]"
    
    print(f"Fetching {family}-tubulin with query: {query}")
    
    raw_entries = fetch_uniprot_json(query, limit=limit)
    entries = [parse_uniprot_entry(e) for e in raw_entries]
    
    # Post-filter by annotation score (can't easily do in query)
    entries = [e for e in entries if e.annotation_score >= min_annotation_score]
    
    print(f"  Got {len(entries)} entries after filtering")
    
    return entries


def write_fasta(entries: list[UniProtEntry], outpath: str):
    with open(outpath, "w") as f:
        for entry in entries:
            f.write(entry.to_fasta() + "\n")


def summarize_entries(entries: list[UniProtEntry]):
    """Print summary statistics."""
    if not entries:
        print("No entries")
        return
    
    reviewed = sum(1 for e in entries if e.reviewed)
    
    # Taxonomic breakdown
    kingdoms = {}
    for e in entries:
        if e.lineage:
            kingdom = e.lineage[0] if e.lineage else "Unknown"
            kingdoms[kingdom] = kingdoms.get(kingdom, 0) + 1
    
    # Annotation score distribution
    scores = {}
    for e in entries:
        scores[e.annotation_score] = scores.get(e.annotation_score, 0) + 1
    
    print(f"\nTotal entries: {len(entries)}")
    print(f"Reviewed (Swiss-Prot): {reviewed}")
    print(f"Unreviewed (TrEMBL): {len(entries) - reviewed}")
    print(f"\nAnnotation score distribution:")
    for score in sorted(scores.keys(), reverse=True):
        print(f"  Score {score}: {scores[score]}")
    print(f"\nTaxonomic distribution (top level):")
    for kingdom, count in sorted(kingdoms.items(), key=lambda x: -x[1]):
        print(f"  {kingdom}: {count}")


# Main execution
if __name__ == "__main__":
    for family in ["alpha", "beta", "gamma", "delta", "epsilon"]:
        print(f"\n{'='*60}")
        print(f"Fetching {family}-tubulin")
        print('='*60)
        
        entries = fetch_tubulin_family(
            family,
            reviewed_only=False,      # Include TrEMBL for coverage
            min_annotation_score=2,   # Relax a bit for rare families
            limit=1000,
        )
        
        summarize_entries(entries)
        
        if entries:
            write_fasta(entries, f"tubulin_{family}_raw.fasta")
            print(f"Wrote tubulin_{family}_raw.fasta")