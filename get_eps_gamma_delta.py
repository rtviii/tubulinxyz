#!/usr/bin/env python3
"""Fetch representative tubulin sequences from UniProt for MSA building."""

import requests
from pathlib import Path

OUTPUT_DIR = Path("/Users/rtviii/dev/tubulinxyz/data/sequences/tubulin")

# Curated UniProt accessions for each family
# Criteria: reviewed (Swiss-Prot), diverse organisms, well-characterized

REPRESENTATIVES = {
    "gamma": [
        # Human and model organisms
        "P23258",   # TUBG1_HUMAN - Human gamma-1
        "P83887",   # TUBG2_HUMAN - Human gamma-2  
        "P53378",   # TBG1_YEAST - S. cerevisiae (TUB4)
        "P33186",   # TBG1_SCHPO - S. pombe
        "P23797",   # TBG_DROME - Drosophila
        "P49553",   # TBG_CAEEL - C. elegans
        "P34787",   # TBG_PLAFO - Plasmodium (parasite)
        "Q55AR3",   # TBG_DICDI - Dictyostelium (amoeba)
        "P27346",   # TBG_ASPNI - Aspergillus (fungus)
    ],
    "delta": [
        # Delta tubulin is rarer - found in organisms with centrioles
        "Q3SYU6",   # TUBD1_HUMAN - Human
        "Q8C0T7",   # TUBD1_MOUSE - Mouse
        "Q5RF10",   # TUBD1_PONAB - Orangutan
        "Q6P9S1",   # TUBD1_DANRE - Zebrafish
        "Q6GQN3",   # TUBD1_XENLA - Xenopus
        "Q9VR85",   # TUBD_DROME - Drosophila
        "Q8I7K3",   # TUBD_PLAF7 - Plasmodium
        "Q23158",   # TUBD_CAEEL - C. elegans
    ],
    "epsilon": [
        # Epsilon tubulin - also centriole-associated
        "Q8N3G8",   # TUBE1_HUMAN - Human
        "Q99MD2",   # TUBE1_MOUSE - Mouse
        "Q6GMB4",   # TUBE1_XENLA - Xenopus
        "Q7ZW40",   # TUBE1_DANRE - Zebrafish  
        "Q9VEN3",   # TUBE_DROME - Drosophila
        "Q8IK81",   # TUBE_PLAF7 - Plasmodium (from your list)
    ],
}


def fetch_uniprot_fasta(accession: str) -> str | None:
    """Fetch a single sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        resp = requests.get(url, timeout=30)
        if resp.status_code == 200:
            return resp.text.strip()
        else:
            print(f"  Warning: {accession} returned status {resp.status_code}")
            return None
    except Exception as e:
        print(f"  Error fetching {accession}: {e}")
        return None


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    for family, accessions in REPRESENTATIVES.items():
        print(f"\n{'='*60}")
        print(f"Fetching {family} tubulin representatives")
        print(f"{'='*60}")
        
        sequences = []
        for acc in accessions:
            print(f"  Fetching {acc}...", end=" ")
            fasta = fetch_uniprot_fasta(acc)
            if fasta:
                sequences.append(fasta)
                # Extract name from header for confirmation
                header = fasta.split("\n")[0]
                print(f"OK - {header[1:60]}...")
            else:
                print("FAILED")
        
        # Write unaligned fasta
        output_file = OUTPUT_DIR / f"tubulin_{family}_representatives.fasta"
        with open(output_file, "w") as f:
            f.write("\n".join(sequences))
        
        print(f"\nWrote {len(sequences)} sequences to {output_file}")


if __name__ == "__main__":
    main()