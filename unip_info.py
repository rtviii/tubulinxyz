#!/usr/bin/env python3
"""Fetch tubulin sequences from UniProt using search queries."""

import requests
from pathlib import Path
import time

OUTPUT_DIR = Path("/Users/rtviii/dev/tubulinxyz/data/sequences/tubulin")

def search_uniprot(query: str, limit: int = 10) -> str:
    """Search UniProt and return FASTA results."""
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "fasta",
        "size": limit,
    }
    
    resp = requests.get(base_url, params=params, timeout=60)
    resp.raise_for_status()
    return resp.text


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Search queries for each family
    # Using reviewed:true to get Swiss-Prot entries only
    queries = {
        "gamma": '(protein_name:"tubulin gamma" OR gene:TUBG1 OR gene:TUBG2 OR gene:TUB4) AND reviewed:true AND length:[400 TO 550]',
        "delta": '(protein_name:"tubulin delta" OR gene:TUBD1) AND reviewed:true AND length:[400 TO 550]',
        "epsilon": '(protein_name:"tubulin epsilon" OR gene:TUBE1) AND reviewed:true AND length:[400 TO 550]',
    }
    
    for family, query in queries.items():
        print(f"\n{'='*60}")
        print(f"Searching for {family} tubulin")
        print(f"Query: {query}")
        print(f"{'='*60}")
        
        try:
            fasta_text = search_uniprot(query, limit=15)
            
            # Count sequences
            n_seqs = fasta_text.count(">")
            print(f"  Found {n_seqs} sequences")
            
            if n_seqs == 0:
                print("  WARNING: No sequences found!")
                continue
            
            # Show headers
            for line in fasta_text.split("\n"):
                if line.startswith(">"):
                    print(f"    {line[:70]}...")
            
            # Save unaligned
            output_file = OUTPUT_DIR / f"tubulin_{family}_representatives.fasta"
            with open(output_file, "w") as f:
                f.write(fasta_text)
            print(f"\n  Saved to {output_file}")
            
        except Exception as e:
            print(f"  ERROR: {e}")
        
        time.sleep(1)  # Be nice to UniProt


if __name__ == "__main__":
    main()