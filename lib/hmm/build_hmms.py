#!/usr/bin/env python3
"""
Build HMMs for all Tubulin and MAP families.
Input:  data/sequences/{type}/{family}.fasta
Output: data/sequences/{type}/{family}.afasta (MSA)
        data/hmms/{type}/{family}.hmm (HMM)
"""

import subprocess
import tempfile
import os
from io import StringIO
from itertools import chain

import pyhmmer
from pyhmmer.easel import Alphabet, TextSequence, TextMSA
from pyhmmer.plan7 import HMM
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger

from lib.models.types_tubulin import TubulinFamily, MapFamily, HmmFamily
from lib.hmm import (
    get_muscle_bin, 
    get_fasta_path, 
    get_aligned_fasta_path, 
    get_hmm_path
)


def muscle_align(seq_records: list[SeqRecord]) -> list[SeqRecord]:
    """Align sequences with MUSCLE."""
    muscle_bin = str(get_muscle_bin())
    
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp:
        SeqIO.write(seq_records, tmp, "fasta")
        tmp_path = tmp.name

    try:
        # Run MUSCLE with quiet flag
        result = subprocess.run(
            [muscle_bin, "-in", tmp_path, "-quiet"],
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            raise RuntimeError(f"MUSCLE failed: {result.stderr}")
        
        aligned = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
        return aligned
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def build_hmm_from_msa(name: str, aligned_seqs: list[SeqRecord], alphabet: Alphabet) -> HMM:
    """Convert aligned SeqRecords to a PyHMMER HMM."""
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


def process_family(family: HmmFamily):
    """Pipeline: Load Fasta -> Align -> Build HMM -> Save Files."""
    fasta_path = get_fasta_path(family)
    afasta_path = get_aligned_fasta_path(family)
    hmm_path = get_hmm_path(family)
    
    if not fasta_path.exists():
        logger.warning(f"[{family.value}] Missing seed FASTA at {fasta_path}")
        return
    
    # 1. Load Sequences
    seqs = list(SeqIO.parse(fasta_path, "fasta"))
    if len(seqs) < 2:
        logger.warning(f"[{family.value}] Too few sequences ({len(seqs)}), skipping.")
        return

    logger.info(f"[{family.value}] Processing {len(seqs)} sequences...")

    # 2. Align (MSA)
    try:
        aligned_seqs = muscle_align(seqs)
        # Save MSA to .afasta
        with open(afasta_path, "w") as f:
            SeqIO.write(aligned_seqs, f, "fasta")
    except Exception as e:
        logger.error(f"[{family.value}] Alignment failed: {e}")
        return

    # 3. Build HMM
    try:
        alphabet = Alphabet.amino()
        hmm = build_hmm_from_msa(family.value, aligned_seqs, alphabet)
        
        # Save HMM
        with open(hmm_path, "wb") as f:
            hmm.write(f)
            
        logger.success(f"[{family.value}] Built HMM (M={hmm.M}, nseq={hmm.nseq}) -> {hmm_path.name}")
    except Exception as e:
        logger.error(f"[{family.value}] HMM build failed: {e}")


def main():
    """Build HMMs for all known families."""
    families = list(chain(TubulinFamily, MapFamily))
    logger.info(f"Starting build for {len(families)} families...")
    
    for family in families:
        process_family(family)

if __name__ == "__main__":
    main()