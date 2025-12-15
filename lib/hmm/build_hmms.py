#!/usr/bin/env python3
"""
Build HMMs for each tubulin family from processed FASTA files.
Usage: python -m lib.hmm.build_hmms
"""

import subprocess
import tempfile
import os
from io import StringIO
from pathlib import Path

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextMSA
from pyhmmer.plan7 import HMM
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from lib.models.types_tubulin import TubulinFamily
from lib.hmm import get_processed_dir, get_hmm_dir, get_muscle_bin, get_seed_fasta_path, get_hmm_path


def muscle_align(seq_records: list[SeqRecord]) -> list[SeqRecord]:
    """Align sequences with MUSCLE, return aligned SeqRecords."""
    muscle_bin = str(get_muscle_bin())
    
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp:
        SeqIO.write(seq_records, tmp, "fasta")
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            [muscle_bin, "-in", tmp_path, "-quiet"],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"MUSCLE failed: {result.stderr}")
        
        aligned = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
        logger.info(f"Aligned {len(aligned)} sequences")
        return aligned
    finally:
        os.unlink(tmp_path)


def build_hmm(name: str, aligned_seqs: list[SeqRecord], alphabet: Alphabet) -> HMM:
    """Build HMM from aligned sequences."""
    text_seqs = [
        TextSequence(name=seq.id.encode(), sequence=str(seq.seq))
        for seq in aligned_seqs
    ]
    
    msa = TextMSA(name.encode(), sequences=text_seqs)
    digital_msa = msa.digitize(alphabet)
    
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    
    hmm, _profile, _optimized = builder.build_msa(digital_msa, background)
    return hmm


def build_family_hmm(family: TubulinFamily) -> HMM | None:
    """Build HMM for a single tubulin family."""
    fasta_path = get_seed_fasta_path(family.value)
    
    if not fasta_path.exists():
        logger.warning(f"Missing FASTA: {fasta_path}")
        return None
    
    logger.info(f"Building HMM for {family.value}...")
    
    seqs = list(SeqIO.parse(fasta_path, "fasta"))
    logger.info(f"  Loaded {len(seqs)} sequences")
    
    if len(seqs) < 2:
        logger.warning(f"  Need at least 2 sequences for MSA, skipping")
        return None
    
    aligned = muscle_align(seqs)
    alphabet = Alphabet.amino()
    hmm = build_hmm(f"tubulin_{family.value}", aligned, alphabet)
    
    logger.info(f"  Built HMM: M={hmm.M}, nseq={hmm.nseq}")
    return hmm


def build_all_hmms():
    """Build HMMs for all tubulin families and save to disk."""
    hmm_dir = get_hmm_dir()
    hmm_dir.mkdir(exist_ok=True)
    
    for family in TubulinFamily:
        hmm = build_family_hmm(family)
        if hmm is None:
            continue
            
        hmm_path = get_hmm_path(family.value)
        with open(hmm_path, "wb") as f:
            hmm.write(f)
        logger.info(f"  Saved to {hmm_path}")


if __name__ == "__main__":
    build_all_hmms()