#!/usr/bin/env python3
"""Diagnose MSA quality for the new clean alignments."""

from pathlib import Path
from collections import Counter
from Bio import AlignIO, SeqIO

MSA_DIR = Path("/Users/rtviii/dev/tubulinxyz/data/sequences/tubulin")

def analyze_msa(msa_path: Path):
    print(f"\n{'='*70}")
    print(f"FILE: {msa_path.name}")
    print(f"{'='*70}")
    
    if not msa_path.exists():
        print("  FILE NOT FOUND")
        return
    
    alignment = AlignIO.read(str(msa_path), "fasta")
    n_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()
    
    print(f"  Sequences: {n_seqs}")
    print(f"  Alignment length: {aln_len}")
    
    # Check sequence IDs
    print(f"\n  Sequences in alignment:")
    for rec in alignment:
        seq_len = len(str(rec.seq).replace("-", ""))
        print(f"    - {rec.id[:50]:<50} ({seq_len} residues)")
    
    # Analyze column quality
    gap_fractions = []
    consensus_strengths = []
    consensus = []
    
    for i in range(aln_len):
        column = alignment[:, i]
        residues = [r for r in column if r not in ("-", ".")]
        
        gap_frac = 1 - (len(residues) / n_seqs)
        gap_fractions.append(gap_frac)
        
        if residues:
            most_common = Counter(residues).most_common(1)[0]
            consensus.append(most_common[0])
            consensus_strength = most_common[1] / len(residues)
        else:
            consensus.append("-")
            consensus_strength = 0
        consensus_strengths.append(consensus_strength)
    
    # Summary stats
    high_gap_cols = sum(1 for g in gap_fractions if g > 0.5)
    very_high_gap_cols = sum(1 for g in gap_fractions if g > 0.9)
    
    avg_consensus = sum(consensus_strengths) / len(consensus_strengths) if consensus_strengths else 0
    weak_consensus_cols = sum(1 for c in consensus_strengths if c < 0.5)
    
    print(f"\n  Gap analysis:")
    print(f"    Columns with >50% gaps: {high_gap_cols} ({100*high_gap_cols/aln_len:.1f}%)")
    print(f"    Columns with >90% gaps: {very_high_gap_cols} ({100*very_high_gap_cols/aln_len:.1f}%)")
    
    print(f"\n  Consensus analysis:")
    print(f"    Average consensus strength: {avg_consensus:.2f}")
    print(f"    Columns with <50% consensus: {weak_consensus_cols} ({100*weak_consensus_cols/aln_len:.1f}%)")
    
    # Show first 80 chars of consensus
    print(f"\n  First 80 consensus positions:")
    print(f"    {''.join(consensus[:80])}")
    
    # Show last 80 chars of consensus  
    print(f"\n  Last 80 consensus positions:")
    print(f"    {''.join(consensus[-80:])}")


def main():
    for family in ["gamma", "delta", "epsilon"]:
        clean_msa = MSA_DIR / f"tubulin_{family}_clean.afasta"
        analyze_msa(clean_msa)
    
    # Also check alpha/beta for comparison
    print("\n\n" + "="*70)
    print("FOR COMPARISON - your working alpha/beta alignments:")
    print("="*70)
    
    for family in ["alpha", "beta"]:
        msa_path = Path(f"/Users/rtviii/dev/tubulinxyz/data/{family}_tubulin/{family}_tubulin.afasta")
        if msa_path.exists():
            alignment = AlignIO.read(str(msa_path), "fasta")
            print(f"\n  {family}: {len(alignment)} seqs, {alignment.get_alignment_length()} columns")


if __name__ == "__main__":
    main()
