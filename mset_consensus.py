import os
from pathlib import Path
from typing import List, Dict
from Bio import AlignIO
import tempfile
import subprocess
import json

from api.config import PROJECT_ROOT
from api.main import MUSCLE_BINARY

class UTNMapper:
    """Maps between MA coordinates and UTN coordinates"""
    
    def __init__(self, master_alignment_path: str, utn_consensus: str, muscle_binary: str):
        self.ma_path = Path(master_alignment_path)
        self.muscle_binary = muscle_binary
        self.utn_consensus = utn_consensus
        
        # Get MA length
        ma_aln = AlignIO.read(str(self.ma_path), "fasta")
        self.ma_length = ma_aln.get_alignment_length()
        
        print(f"UTN consensus length: {len(self.utn_consensus)}")
        print(f"MA length: {self.ma_length}")
        
        # Initialize mappings
        self.ma_to_utn: List[int] = []  # Index = MA position (0-based), Value = UTN position (1-based) or -1
        self.utn_to_ma: List[int] = []  # Index = UTN position (0-based), Value = MA position (1-based) or -1/-2
        
        # Perform alignment and build mappings
        self._align_and_map()
    
    def _align_and_map(self):
        """Profile align UTN consensus to MA and build bidirectional mappings"""
        
        # Write UTN consensus to temp file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(">UTN_CONSENSUS\n")
            f.write(self.utn_consensus + "\n")
            utn_temp = f.name
        
        # Run MUSCLE profile alignment
        with tempfile.NamedTemporaryFile(mode='w', suffix='.aln', delete=False) as f:
            output_temp = f.name
        
        try:
            cmd = [
                self.muscle_binary, '-profile',
                '-in1', str(self.ma_path),
                '-in2', utn_temp,
                '-out', output_temp
            ]
            
            print("Running MUSCLE profile alignment...")
            subprocess.run(cmd, check=True, capture_output=True)
            
            # Parse the alignment
            self._parse_alignment(output_temp)
            
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"MUSCLE Error: {e.stderr.decode()}")
        finally:
            Path(utn_temp).unlink(missing_ok=True)
            Path(output_temp).unlink(missing_ok=True)
    
    def _parse_alignment(self, alignment_file: str):
        """Parse the profile alignment and build mappings"""
        
        alignment = AlignIO.read(alignment_file, "fasta")
        
        # Find the UTN consensus in the alignment
        utn_record = next((r for r in alignment if r.id == "UTN_CONSENSUS"), None)
        if not utn_record:
            raise ValueError("UTN_CONSENSUS not found in alignment output")
        
        aligned_utn = str(utn_record.seq)
        
        # Get MA sequences (all except UTN_CONSENSUS)
        ma_sequences = [str(r.seq) for r in alignment if r.id != "UTN_CONSENSUS"]
        
        print(f"Aligned length: {len(aligned_utn)}")
        
        # Determine column types: 'original' MA column vs 'insertion' column
        # Original = at least one MA sequence has a residue
        # Insertion = ALL MA sequences have gaps (new column added by MUSCLE)
        ma_column_types = []
        for i in range(len(aligned_utn)):
            has_ma_residue = any(seq[i] != '-' for seq in ma_sequences)
            ma_column_types.append('original' if has_ma_residue else 'insertion')
        
        # Build mappings
        ma_index = 0   # Counter for original MA positions (0-based)
        utn_index = 0  # Counter for UTN positions (0-based)
        
        for col_idx, (utn_char, col_type) in enumerate(zip(aligned_utn, ma_column_types)):
            
            if col_type == 'insertion':
                # New column: UTN has insertion relative to MA
                if utn_char != '-':
                    self.utn_to_ma.append(-2)  # -2 = insertion
                    utn_index += 1
                    
            else:  # col_type == 'original'
                # Original MA column
                current_ma_pos = ma_index + 1  # 1-based
                
                if utn_char == '-':
                    # MA has this position, UTN doesn't
                    self.ma_to_utn.append(-1)
                else:
                    # Both have residues (including 'x')
                    current_utn_pos = utn_index + 1  # 1-based
                    self.ma_to_utn.append(current_utn_pos)
                    self.utn_to_ma.append(current_ma_pos)
                    utn_index += 1
                
                ma_index += 1
        
        print(f"MA positions: {len(self.ma_to_utn)}")
        print(f"UTN positions: {len(self.utn_to_ma)}")
        print(f"UTN insertions (no MA equivalent): {self.utn_to_ma.count(-2)}")
        print(f"MA missing (no UTN equivalent): {self.ma_to_utn.count(-1)}")
    
    def ma_to_utn_pos(self, ma_position: int) -> int:
        """Convert MA position (1-based) to UTN position (1-based), or -1 if no mapping"""
        if 1 <= ma_position <= len(self.ma_to_utn):
            return self.ma_to_utn[ma_position - 1]
        return -1
    
    def utn_to_ma_pos(self, utn_position: int) -> int:
        """Convert UTN position (1-based) to MA position (1-based), -1 if no mapping, -2 if insertion"""
        if 1 <= utn_position <= len(self.utn_to_ma):
            return self.utn_to_ma[utn_position - 1]
        return -1
    
    def save_mappings(self, output_path: str):
        """Save mappings to JSON"""
        output = {
            "ma_to_utn": self.ma_to_utn,
            "utn_to_ma": self.utn_to_ma,
            "stats": {
                "ma_length": len(self.ma_to_utn),
                "utn_length": len(self.utn_to_ma),
                "utn_insertions": self.utn_to_ma.count(-2),
                "ma_missing": self.ma_to_utn.count(-1),
                "mapped_positions": len([x for x in self.ma_to_utn if x != -1])
            }
        }
        
        with open(output_path, "w") as f:
            json.dump(output, f, indent=2)
        
        print(f"\nMappings saved to {output_path}")


def main():
    # UTN Consensus from the paper
    UTN_CONSENSUS = "MRECISIHIGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGKHVPRAVFVDLEPTVIDEVRTGTYRQLFHPEQLISGKEDAANNYARGHYTIGKEIVDLxLDRIRKLADNCTGLQGFLVFHSxGGGTGSGxGSLLMERLSVDYGKKSKLxFTIYPSPQVSTAVVEPYNSVLTTHTxLEHTDxAxMVDNEAIYDICRRNLDIERPTYTNLNRLISQVISSLTASLRFDGALNVDLTEFQTNLVPYPRIHFxLSSYAPVISAEKAYHEQLSVAEITNAcFEPANxMVKCDPRHGKYMACCLMYRGDVVPKDVNAAVATIKTKRTIQFVDWCPTGFKxGINYQPPTVVPGGDLAKVQRAVCMLSNTTAIAEaWSRLDHKFDLMYAKRAFVHWYVGEGMEEGEFSEAREDLAALEKDYEEVGADSx"
    
    # Configuration
    MA_PATH = os.path.join(PROJECT_ROOT, 'data/alpha_tubulin/alpha_tubulin.afasta')
    
    # Create mapper
    mapper = UTNMapper(MA_PATH, UTN_CONSENSUS, MUSCLE_BINARY)
    
    # Example queries
    print("\n" + "="*60)
    print("Example mappings:")
    print("="*60)
    for test_pos in [1, 40, 100, 200, 300, 400]:
        utn_pos = mapper.ma_to_utn_pos(test_pos)
        print(f"MA position {test_pos:3d} -> UTN position {utn_pos:3d}")
    
    print()
    for test_pos in [1, 40, 100, 200, 300, 400]:
        ma_pos = mapper.utn_to_ma_pos(test_pos)
        status = "INSERTION" if ma_pos == -2 else "MISSING" if ma_pos == -1 else f"MA {ma_pos}"
        print(f"UTN position {test_pos:3d} -> {status}")
    
    # Save mappings
    mapper.save_mappings("utn_alpha_mapping.json")


if __name__ == "__main__":
    main()