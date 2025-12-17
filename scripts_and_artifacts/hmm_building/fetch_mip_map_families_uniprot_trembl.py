"""
MAP/MIP Sequence Fetcher - Simplified
Just uses gene names to query UniProt. No InterPro complexity.
"""
import requests
import time
import os
from dataclasses import dataclass


MAP_FAMILIES = {
    # ===================
    # Classical Structural MAPs
    # ===================
    "classical_map_tau": {
        "genes": ["MAPT", "MAP2", "MAP4"],
        "min_len": 200, "max_len": 3000,
        "desc": "Tau/MAP2/MAP4 family"
    },
    "map1_family": {
        "genes": ["MAP1A", "MAP1B", "MAP1S", "MAP1LC3A", "MAP1LC3B"],
        "min_len": 100, "max_len": 3500,
        "desc": "MAP1 family"
    },
    "map7_ensconsin": {
        "genes": ["MAP7", "MAP7D1", "MAP7D2", "MAP7D3", "ENSCONSIN"],
        "min_len": 400, "max_len": 1200,
        "desc": "MAP7/Ensconsin family"
    },
    "map6_stop": {
        "genes": ["MAP6", "MAP6D1", "STOP"],
        "min_len": 400, "max_len": 1200,
        "desc": "MAP6/STOP cold-stable MAPs"
    },
    # "map9": {
    #     "genes": ["MAP9", "ASAP"],
    #     "min_len": 400, "max_len": 900,
    #     "desc": "MAP9/ASAP"
    # },
    "doublecortin": {
        "genes": ["DCX", "DCLK1", "DCLK2", "DCLK3"],
        "min_len": 300, "max_len": 1500,
        "desc": "Doublecortin family"
    },

    # ===================
    # Bundlers/Crosslinkers
    # ===================
    "bundler_prc1": {
        "genes": ["PRC1", "ASE1", "MAP65"],
        "min_len": 400, "max_len": 1200,
        "desc": "PRC1/MAP65 crosslinkers"
    },
    "spectraplakin": {
        "genes": ["MACF1", "DST", "BPAG1", "ACF7"],
        "min_len": 2000, "max_len": 8000,
        "desc": "Spectraplakins (MT-actin crosslinkers)"
    },

    # ===================
    # Plus-End Tracking (+TIPs)
    # ===================
    "plus_eb_family": {
        "genes": ["MAPRE1", "MAPRE2", "MAPRE3", "EB1", "BIM1"],
        "min_len": 150, "max_len": 600,
        "desc": "EB1/EB3 family"
    },
    "plus_xmap215_tog": {
        "genes": ["CKAP5", "XMAP215", "CLASP1", "CLASP2", "STU2", "chTOG"],
        "min_len": 800, "max_len": 5000,
        "desc": "TOG domain polymerases"
    },
    "plus_clip": {
        "genes": ["CLIP1", "CLIP2", "CLIP170"],
        "min_len": 300, "max_len": 2000,
        "desc": "CLIP-170 family"
    },
    "plus_tacc": {
        "genes": ["TACC1", "TACC2", "TACC3"],
        "min_len": 400, "max_len": 1000,
        "desc": "TACC family (chTOG partners)"
    },

    # ===================
    # Minus-End
    # ===================
    "minus_camsap": {
        "genes": ["CAMSAP1", "CAMSAP2", "CAMSAP3", "PTRN1", "PATRONIN"],
        "min_len": 400, "max_len": 2500,
        "desc": "Patronin/CAMSAP family"
    },
    "minus_gamma_turc": {
        "genes": ["TUBGCP2", "TUBGCP3", "TUBGCP4", "TUBGCP5", "TUBGCP6"],
        "min_len": 400, "max_len": 2000,
        "desc": "Gamma-tubulin complex proteins"
    },

    # ===================
    # Spindle/Mitotic MAPs
    # ===================
    "spindle_tpx2": {
        "genes": ["TPX2"],
        "min_len": 400, "max_len": 900,
        "desc": "TPX2 spindle assembly factor"
    },
    "spindle_numa": {
        "genes": ["NUMA1", "NuMA", "LIN5"],
        "min_len": 1500, "max_len": 2500,
        "desc": "NuMA spindle pole protein"
    },
    "spindle_astrin": {
        "genes": ["SPAG5", "ASTRIN", "KNSTRN", "SKAP"],
        "min_len": 200, "max_len": 1400,
        "desc": "Astrin-SKAP complex"
    },

    # ===================
    # Destabilizers
    # ===================
    "destab_stathmin": {
        "genes": ["STMN1", "STMN2", "STMN3", "STMN4"],
        "min_len": 100, "max_len": 400,
        "desc": "Stathmin family"
    },
    "destab_severing": {
        "genes": ["KATNA1", "KATNB1", "SPAST", "FIGN", "KATNAL1", "KATNAL2"],
        "min_len": 300, "max_len": 1500,
        "desc": "Severing enzymes (Katanin/Spastin/Fidgetin)"
    },
    "destab_kinesin13": {
        "genes": ["KIF2A", "KIF2B", "KIF2C", "MCAK"],
        "min_len": 400, "max_len": 1200,
        "desc": "Kinesin-13 depolymerases"
    },

    # ===================
    # Tubulin Code - Writers
    # ===================
    "code_ttll": {
        "genes": ["TTLL1", "TTLL2", "TTLL4", "TTLL5", "TTLL6", "TTLL7", "TTLL9", "TTLL11", "TTLL12", "TTLL13"],
        "min_len": 300, "max_len": 1500,
        "desc": "TTLL polyglutamylases/glycylases"
    },
    "code_atat": {
        "genes": ["ATAT1", "MEC17"],
        "min_len": 150, "max_len": 600,
        "desc": "Alpha-tubulin acetyltransferase"
    },

    # ===================
    # Tubulin Code - Erasers
    # ===================
    "code_ccp_deglutamylase": {
        "genes": ["AGTPBP1", "CCP1", "CCP2", "CCP3", "CCP4", "CCP5", "CCP6", "AGBL1", "AGBL2", "AGBL3", "AGBL4", "AGBL5"],
        "min_len": 400, "max_len": 1400,
        "desc": "CCP deglutamylases"
    },
    "code_vash_detyrosinase": {
        "genes": ["VASH1", "VASH2", "SVBP"],
        "min_len": 100, "max_len": 400,
        "desc": "Vasohibin detyrosinases + SVBP"
    },
    "code_hdac6_deacetylase": {
        "genes": ["HDAC6"],
        "min_len": 800, "max_len": 1400,
        "desc": "HDAC6 tubulin deacetylase"
    },
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
            
            # Check for next page
            link = resp.headers.get("Link", "")
            if 'rel="next"' in link:
                next_url = link.split(";")[0].strip("<>")
                base_url = next_url
                params = {}
                time.sleep(0.3)
            else:
                break
        except requests.exceptions.RequestException as e:
            print(f"  Error: {e}")
            break
    
    return results[:limit]


def parse_entry(raw: dict) -> UniProtEntry:
    """Parse UniProt JSON into our dataclass."""
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


def fetch_family(key: str, config: dict, limit: int = 200) -> list[UniProtEntry]:
    """Fetch sequences for a MAP family."""
    genes = config["genes"]
    min_len, max_len = config["min_len"], config["max_len"]
    
    gene_query = " OR ".join(f"gene:{g}" for g in genes)
    query = f"({gene_query}) AND length:[{min_len} TO {max_len}]"
    
    print(f"Fetching {key} ({config['desc']})...")
    
    raw = fetch_uniprot(query, limit=limit)
    entries = [parse_entry(e) for e in raw]
    
    return entries


def write_fasta(entries: list[UniProtEntry], path: str):
    with open(path, "w") as f:
        for e in entries:
            f.write(e.to_fasta() + "\n")


if __name__ == "__main__":
    os.makedirs("map_sequences", exist_ok=True)
    
    print(f"{'='*50}")
    print(f"MAP/MIP Sequence Fetch - {len(MAP_FAMILIES)} families")
    print(f"{'='*50}\n")
    
    for key, config in MAP_FAMILIES.items():
        entries = fetch_family(key, config, limit=1000)
        
        if entries:
            out_file = f"map_sequences/map_{key}.fasta"
            write_fasta(entries, out_file)
            
            reviewed = sum(1 for e in entries if e.reviewed)
            print(f"  -> {len(entries)} sequences ({reviewed} reviewed)")
        else:
            print(f"  -> WARNING: No results")
        
        time.sleep(0.5)
    
    print("\nDone. Check 'map_sequences/' folder.")

