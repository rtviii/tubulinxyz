#!/usr/bin/env python3
"""Clean up the fetched sequences and align them."""

from pathlib import Path
from Bio import SeqIO
import subprocess

DATA_DIR = Path("/Users/rtviii/dev/tubulinxyz/data/sequences/tubulin")
MUSCLE = Path("/Users/rtviii/dev/tubulinxyz/muscle3.8.1")

# Patterns to EXCLUDE from each family
EXCLUDE_PATTERNS = {
    "gamma": ["TBB", "TBA"],  # Exclude beta and alpha tubulins
    "delta": [],              # All look good
    "epsilon": ["TEDC"],      # Exclude complex proteins, keep only TBE
}

def clean_and_align(family: str):
    input_file = DATA_DIR / f"tubulin_{family}_representatives.fasta"
    cleaned_file = DATA_DIR / f"tubulin_{family}_cleaned.fasta"
    aligned_file = DATA_DIR / f"tubulin_{family}_clean.afasta"
    
    print(f"\n{'='*60}")
    print(f"Processing {family} tubulin")
    print(f"{'='*60}")
    
    # Read and filter
    exclude = EXCLUDE_PATTERNS.get(family, [])
    kept = []
    removed = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        # Check if any exclude pattern is in the ID
        should_exclude = any(pat in record.id for pat in exclude)
        
        if should_exclude:
            removed.append(record.id)
        else:
            kept.append(record)
    
    print(f"  Kept: {len(kept)} sequences")
    for r in kept:
        print(f"    + {r.id[:50]}")
    
    if removed:
        print(f"  Removed: {len(removed)} sequences")
        for rid in removed:
            print(f"    - {rid[:50]}")
    
    # Write cleaned fasta
    SeqIO.write(kept, cleaned_file, "fasta")
    
    # Align with MUSCLE
    print(f"\n  Aligning with MUSCLE...")
    cmd = [str(MUSCLE), "-in", str(cleaned_file), "-out", str(aligned_file)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr[:200]}")
        return
    
    # Quick stats on result
    from Bio import AlignIO
    alignment = AlignIO.read(str(aligned_file), "fasta")
    print(f"  Result: {len(alignment)} sequences, {alignment.get_alignment_length()} columns")
    
    print(f"  Saved to: {aligned_file}")


def main():
    for family in ["gamma", "delta", "epsilon"]:
        clean_and_align(family)


if __name__ == "__main__":
    main()