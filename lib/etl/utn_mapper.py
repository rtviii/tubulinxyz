# lib/etl/utn_mapper.py
"""
Maps between Morisette Universal Tubulin Numbering (UTN) coordinates and our
master alignment (MA) indices. Created specifically for Morisette database
ingestion (Morisette et al., PLOS ONE 2023, doi:10.1371/journal.pone.0295279).

The UTN system uses gapless consensus sequences (440 residues each, truncated
at the C-terminal tail) derived from Clustal Omega alignments of ~90 alpha/beta
and ~21 gamma tubulin sequences. UTN position = 1-based index into the consensus.

This module profile-aligns each UTN consensus against our family MSA using MUSCLE,
producing a bidirectional utn_position <-> master_index mapping.

Usage (CLI):
    python -m lib.etl.utn_mapper --family alpha --output mappings.json
    python -m lib.etl.utn_mapper --family beta --output mappings.json
"""

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Optional, List

from Bio import AlignIO
from loguru import logger

from api.config import settings


# ---------------------------------------------------------------------------
# UTN consensus sequences from Morisette et al. S1 Table (verified 2026-03-30)
# Uppercase = strong consensus, lowercase = weak consensus, x = no consensus
# All are exactly 440 characters (CTT truncated)
# ---------------------------------------------------------------------------

UTN_ALPHA_CONSENSUS = (
    "MRECISIHIGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGKHVPRAVFVDLEPTVIDEVRT"
    "GTYRQLFHPEQLISGKEDAANNYARGHYTIGKEIVDLxLDRIRKLADNCTGLQGFLVFHSxGGGTGSGxGSLLMERLSVD"
    "YGKKSKLxFTIYPSPQVSTAVVEPYNSVLTTHTxLEHTDxAxMVDNEAIYDICRRNLDIERPTYTNLNRLISQVISSLTA"
    "SLRFDGALNVDLTEFQTNLVPYPRIHFxLSSYAPVISAEKAYHEQLSVAEITNAcFEPANxMVKCDPRHGKYMACCLMYR"
    "GDVVPKDVNAAVATIKTKRTIQFVDWCPTGFKxGINYQPPTVVPGGDLAKVQRAVCMLSNTTAIAEaWSRLDHKFDLMYA"
    "KRAFVHWYVGEGMEEGEFSEAREDLAALEKDYEEVGADSx"
)

UTN_BETA_CONSENSUS = (
    "MREIVHIQxGQCGNQIGAKFWEVIxDEHGIDPTGxYxGDSDLQLERINVYYNEASGGRYVPRAVLMDLEPGTMDSVRSGP"
    "FGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELIDSVLDVVRKEAENCDCLQGFQIxHSLGGGTGSGMGTLLISKIREEYP"
    "DRMMxTFSVxPSPKVSDTVVEPYNATLSVHQLVENADExfxIDNEALYDICFRTLKLTTPTYGDLNHLVSATMSGVTxCL"
    "RFPGQLNSDLRKLAVNMVPFPRLHFFMxGFAPLTSRGSQQYRALTVPELTQQMWDAKNMMcAxDPRHGRYLTAAAMFRGR"
    "MSTKEVDEQMLNVQNKNSSYFVEWIPNNVKTSVCDIPPRGLKMAATFIGNSTAIQEMFKRVSEQFTAMFRRKAFLHWYTG"
    "EGMDEMEFTEAESNMNDLVSEYQQYQDATADxEEEExEEE"
)

UTN_GAMMA_CONSENSUS = (
    "MPREIITLQVGQCGNQIGxEFWKQLCxEHGISPEGILEDFATEGxDRKDVFFYQADDEHYIPRAILIDLEPRVINxIQNS"
    "xYSxLYNPENIFISKHGGGAGNNWASGYSQGEKVQEDIxDMIDREADGSDSLEGFVLCHSIAGGTGSGMGSYLLERLNDR"
    "YPKKLIQTYSVFPNQxExSDVVVQPYNSLLTLKRLTQNADCVVVLDNTALNRIAxDRLHIxNPTFSQxNSLVSTVMSAST"
    "TTLRYPGYMNNDLVGLIASLIPTPRCHFLMTGYTPLTxDQxVSSVRKTTVLDVMRRLLQPKNIMVSTxxKxxxNxxYISI"
    "LNIIQGEVDPTQVHKSLQRIRERKLANFIPWGPASIQVALSRKSPYIQTSHRVSGLMLANHTSISSLFERxLxQYDKLRK"
    "RxAFLDQYRKExMFKDNLDEFDESREVVQxLIDEYKAAER"
)

# Family -> (consensus, MSA filename relative to data/)
_FAMILY_CONFIG = {
    "tubulin_alpha": (UTN_ALPHA_CONSENSUS, "alpha_tubulin/alpha_tubulin.afasta"),
    "tubulin_beta":  (UTN_BETA_CONSENSUS,  "beta_tubulin/beta_tubulin.afasta"),
    "tubulin_gamma": (UTN_GAMMA_CONSENSUS, "gamma_tubulin/tubulin_gamma_clean.afasta"),
}


class UTNMapper:
    """Maps between Morisette UTN coordinates and our master alignment indices.

    Profile-aligns the UTN consensus sequence against our family MSA using MUSCLE,
    producing a bidirectional utn_position <-> master_index mapping.
    """

    def __init__(
        self,
        family: str,
        utn_consensus: str,
        msa_path: Path,
        muscle_binary: str,
    ):
        self.family = family
        self.utn_consensus = utn_consensus
        self.msa_path = Path(msa_path)
        self.muscle_binary = muscle_binary

        # utn_to_master_map[utn_pos] = master_index or None
        # master_to_utn_map[master_idx] = utn_pos or None
        self.utn_to_master_map: Dict[int, Optional[int]] = {}
        self.master_to_utn_map: Dict[int, Optional[int]] = {}

        self._stats: Dict[str, int] = {}

        self._align_and_map()

    def _align_and_map(self):
        """Profile-align UTN consensus to our MSA and build bidirectional mappings."""

        seq_id = "UTN_CONSENSUS"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(f">{seq_id}\n{self.utn_consensus}\n")
            seq_temp = Path(f.name)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".aln", delete=False) as f:
            out_temp = Path(f.name)

        try:
            cmd = [
                self.muscle_binary, "-profile",
                "-in1", str(self.msa_path),
                "-in2", str(seq_temp),
                "-out", str(out_temp),
            ]
            logger.info(f"Running MUSCLE profile alignment for {self.family}...")
            subprocess.run(cmd, check=True, capture_output=True)

            alignment = AlignIO.read(str(out_temp), "fasta")
            utn_record = next((r for r in alignment if r.id == seq_id), None)
            if not utn_record:
                raise ValueError("UTN_CONSENSUS not found in alignment output")

            aligned_utn = str(utn_record.seq)
            ma_sequences = [str(r.seq) for r in alignment if r.id != seq_id]
            n_cols = len(aligned_utn)

            # Classify each column: "original" (exists in our MSA) or "insertion"
            # (added by MUSCLE to accommodate UTN insertions relative to our MSA)
            is_original = []
            for i in range(n_cols):
                has_ma_residue = any(seq[i] not in ("-", ".") for seq in ma_sequences)
                is_original.append(has_ma_residue)

            # Walk columns and build mappings
            ma_idx = 0   # 0-based counter for original MA columns
            utn_idx = 0  # 0-based counter for UTN residues

            utn_insertions = 0
            ma_gaps = 0
            mapped = 0

            for col_i in range(n_cols):
                utn_char = aligned_utn[col_i]

                if not is_original[col_i]:
                    # Insertion column (not in our MSA)
                    if utn_char not in ("-", "."):
                        utn_pos = utn_idx + 1
                        self.utn_to_master_map[utn_pos] = None
                        utn_insertions += 1
                        utn_idx += 1
                else:
                    # Original MA column
                    master_pos = ma_idx + 1

                    if utn_char in ("-", "."):
                        # MA has this position, UTN doesn't
                        self.master_to_utn_map[master_pos] = None
                        ma_gaps += 1
                    else:
                        # Both have residues -- mapped
                        utn_pos = utn_idx + 1
                        self.utn_to_master_map[utn_pos] = master_pos
                        self.master_to_utn_map[master_pos] = utn_pos
                        mapped += 1
                        utn_idx += 1

                    ma_idx += 1

            self._stats = {
                "family": self.family,
                "utn_length": len(self.utn_consensus),
                "ma_length": ma_idx,
                "mapped": mapped,
                "utn_insertions": utn_insertions,
                "ma_gaps": ma_gaps,
            }

            logger.info(
                f"UTNMapper [{self.family}]: {mapped} mapped, "
                f"{utn_insertions} UTN insertions, {ma_gaps} MA-only positions"
            )

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MUSCLE failed: {e.stderr.decode()}")
        finally:
            seq_temp.unlink(missing_ok=True)
            out_temp.unlink(missing_ok=True)

    def utn_to_master(self, utn_pos: int) -> Optional[int]:
        """Convert UTN position (1-based) to master alignment index, or None if unmapped."""
        return self.utn_to_master_map.get(utn_pos)

    def master_to_utn(self, master_pos: int) -> Optional[int]:
        """Convert master alignment index (1-based) to UTN position, or None if unmapped."""
        return self.master_to_utn_map.get(master_pos)

    @property
    def stats(self) -> Dict[str, int]:
        return self._stats

    def save_mappings(self, output_path: str):
        """Save mappings to JSON for inspection/debugging."""
        output = {
            "stats": self._stats,
            "utn_to_master": {str(k): v for k, v in sorted(self.utn_to_master_map.items())},
            "master_to_utn": {str(k): v for k, v in sorted(self.master_to_utn_map.items())},
        }
        with open(output_path, "w") as f:
            json.dump(output, f, indent=2)
        logger.info(f"Mappings saved to {output_path}")


def get_mapper(family: str, muscle_binary: str = None) -> UTNMapper:
    """Factory: create a UTNMapper for a given family using default paths."""
    if family not in _FAMILY_CONFIG:
        raise ValueError(f"Unknown family '{family}'. Available: {list(_FAMILY_CONFIG.keys())}")

    consensus, msa_rel_path = _FAMILY_CONFIG[family]
    msa_path = settings.PROJECT_ROOT / "data" / msa_rel_path

    if not msa_path.exists():
        raise FileNotFoundError(f"MSA file not found: {msa_path}")

    if muscle_binary is None:
        muscle_binary = settings.MUSCLE_BINARY

    return UTNMapper(family, consensus, msa_path, muscle_binary)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate UTN <-> MA mappings for a tubulin family")
    parser.add_argument("--family", required=True, choices=list(_FAMILY_CONFIG.keys()),
                        help="Tubulin family to map")
    parser.add_argument("--output", default=None,
                        help="Output JSON path (default: prints stats only)")
    parser.add_argument("--muscle", default=None,
                        help="Path to MUSCLE binary (default: from settings)")
    args = parser.parse_args()

    mapper = get_mapper(args.family, muscle_binary=args.muscle)

    print(f"\nStats: {json.dumps(mapper.stats, indent=2)}")

    # Spot-check a few positions
    print("\nSample mappings (UTN -> MA):")
    for pos in [1, 40, 100, 200, 300, 400, 440]:
        ma = mapper.utn_to_master(pos)
        print(f"  UTN {pos:3d} -> MA {ma}")

    if args.output:
        mapper.save_mappings(args.output)
