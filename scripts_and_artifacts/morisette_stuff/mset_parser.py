import os
import pandas as pd
import json
import re
from typing import Dict, List, Any

from api.config import PROJECT_ROOT
from api.main import MUSCLE_BINARY

class TubulinDataParser:
    
    def __init__(self, utn_mapper):
        """
        utn_mapper: UTNMapper instance with ma_to_utn and utn_to_ma mappings
        """
        self.mapper = utn_mapper
    
    def parse_mutation_aa1(self, aa1: str) -> Dict[str, Any]:
        """
        Parse mutation notation like 'M1L' or 'R2A'
        Returns: {wt: 'M', position: 1, mut: 'L'}
        """
        match = re.match(r'([A-Z])(\d+)([A-Z])', aa1)
        if match:
            return {
                'wild_type': match.group(1),
                'utn_position': int(match.group(2)),
                'mutant': match.group(3)
            }
        return None
    
    def parse_modification_aa1(self, aa1: str) -> Dict[str, Any]:
        """
        Parse modification notation like 'C4-PLM' or 'S54-PHO'
        Returns: {aa: 'C', position: 4, mod_type: 'PLM'}
        """
        match = re.match(r'([A-Z])(\d+)-([A-Z]+)', aa1)
        if match:
            return {
                'amino_acid': match.group(1),
                'utn_position': int(match.group(2)),
                'modification_type': match.group(3)
            }
        return None
    
    def parse_mutations(self, csv_path: str, output_json: str = None) -> List[Dict]:
        """
        Parse alpha-tubulin mutations CSV
        """
        # Read CSV, skip the first 2 header rows
        df = pd.read_csv(csv_path, skiprows=1)
        
        # Filter for alpha tubulins only
        df = df[df['UTN'].str.contains('α', na=False)]
        
        mutations = []
        
        for _, row in df.iterrows():
            aa1 = str(row['AA1']).strip()
            mutation_info = self.parse_mutation_aa1(aa1)
            
            if not mutation_info:
                print(f"Warning: Could not parse AA1: {aa1}")
                continue
            
            utn_pos = mutation_info['utn_position']
            ma_pos = self.mapper.utn_to_ma_pos(utn_pos)
            
            entry = {
                'utn_position': utn_pos,
                'ma_position': ma_pos,
                'ma_status': 'mapped' if ma_pos > 0 else ('insertion' if ma_pos == -2 else 'missing'),
                'wild_type': mutation_info['wild_type'],
                'mutant': mutation_info['mutant'],
                'tubulin_type': str(row['UTN']).strip(),
                'tubulin_id': str(row['Tubulin']).strip(),
                'species': str(row['Species']).strip(),
                'phenotype': str(row['Phenotype']).strip(),
                'reference': str(row['Reference']).strip() if 'Reference' in row else '',
                'ref_link': str(row['Ref Link Info']).strip() if 'Ref Link Info' in row else '',
                'keywords': str(row['Keywords']).strip() if 'Keywords' in row else '',
                'notes': str(row['Notes']).strip() if 'Notes' in row else ''
            }
            
            mutations.append(entry)
        
        print(f"Parsed {len(mutations)} alpha-tubulin mutations")
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(mutations, f, indent=2)
            print(f"Saved to {output_json}")
        
        return mutations
    
    def parse_modifications(self, csv_path: str, output_json: str = None) -> List[Dict]:
        """
        Parse alpha-tubulin modifications CSV
        """
        # Read CSV, skip the first 2 header rows
        df = pd.read_csv(csv_path, skiprows=1)
        
        # Filter for alpha tubulins only
        df = df[df['UTN'].str.contains('α', na=False)]
        
        modifications = []
        
        for _, row in df.iterrows():
            aa1 = str(row['AA1']).strip()
            mod_info = self.parse_modification_aa1(aa1)
            
            if not mod_info:
                print(f"Warning: Could not parse AA1: {aa1}")
                continue
            
            utn_pos = mod_info['utn_position']
            ma_pos = self.mapper.utn_to_ma_pos(utn_pos)
            
            entry = {
                'utn_position': utn_pos,
                'ma_position': ma_pos,
                'ma_status': 'mapped' if ma_pos > 0 else ('insertion' if ma_pos == -2 else 'missing'),
                'amino_acid': mod_info['amino_acid'],
                'modification_type': mod_info['modification_type'],
                'tubulin_type': str(row['UTN']).strip(),
                'tubulin_id': str(row['Tubulin']).strip(),
                'species': str(row['Species']).strip(),
                'phenotype': str(row['Phenotype']).strip(),
                'database': str(row['Database']).strip() if 'Database' in row else '',
                'db_link': str(row['DBLink:']).strip() if 'DBLink:' in row else '',
                'keywords': str(row['Key words:']).strip() if 'Key words:' in row else '',
                'notes': str(row['Notes:']).strip() if 'Notes:' in row else ''
            }
            
            modifications.append(entry)
        
        print(f"Parsed {len(modifications)} alpha-tubulin modifications")
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(modifications, f, indent=2)
            print(f"Saved to {output_json}")
        
        return modifications
    
    def get_statistics(self, mutations: List[Dict], modifications: List[Dict]) -> Dict:
        """
        Get summary statistics
        """
        mut_mapped = len([m for m in mutations if m['ma_status'] == 'mapped'])
        mut_insertions = len([m for m in mutations if m['ma_status'] == 'insertion'])
        mut_missing = len([m for m in mutations if m['ma_status'] == 'missing'])
        
        mod_mapped = len([m for m in modifications if m['ma_status'] == 'mapped'])
        mod_insertions = len([m for m in modifications if m['ma_status'] == 'insertion'])
        mod_missing = len([m for m in modifications if m['ma_status'] == 'missing'])
        
        return {
            'mutations': {
                'total': len(mutations),
                'mapped_to_ma': mut_mapped,
                'utn_insertions': mut_insertions,
                'missing_in_utn': mut_missing
            },
            'modifications': {
                'total': len(modifications),
                'mapped_to_ma': mod_mapped,
                'utn_insertions': mod_insertions,
                'missing_in_utn': mod_missing
            }
        }


def main():
    from mset_consensus import UTNMapper
    
    # Configuration
    UTN_CONSENSUS = "MRECISIHIGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGKHVPRAVFVDLEPTVIDEVRTGTYRQLFHPEQLISGKEDAANNYARGHYTIGKEIVDLxLDRIRKLADNCTGLQGFLVFHSxGGGTGSGxGSLLMERLSVDYGKKSKLxFTIYPSPQVSTAVVEPYNSVLTTHTxLEHTDxAxMVDNEAIYDICRRNLDIERPTYTNLNRLISQVISSLTASLRFDGALNVDLTEFQTNLVPYPRIHFxLSSYAPVISAEKAYHEQLSVAEITNAcFEPANxMVKCDPRHGKYMACCLMYRGDVVPKDVNAAVATIKTKRTIQFVDWCPTGFKxGINYQPPTVVPGGDLAKVQRAVCMLSNTTAIAEaWSRLDHKFDLMYAKRAFVHWYVGEGMEEGEFSEAREDLAALEKDYEEVGADSx"
    MA_PATH = os.path.join(PROJECT_ROOT, 'data/alpha_tubulin/alpha_tubulin.afasta')
    mapper = UTNMapper(MA_PATH, UTN_CONSENSUS, MUSCLE_BINARY)
    
    MUTATIONS_CSV     = "data/alpha_tubulin/alpha_tubulin_mutations.csv"
    MODIFICATIONS_CSV = "data/alpha_tubulin/alpha_tubulin_modifications.csv"
    
    # Create mapper
    print("Creating UTN to MA mapper...")
    mapper = UTNMapper(MA_PATH, UTN_CONSENSUS, MUSCLE_BINARY)
    
    # Parse data
    parser = TubulinDataParser(mapper)
    
    print("\nParsing mutations...")
    mutations = parser.parse_mutations(MUTATIONS_CSV, "alpha_mutations.json")
    
    print("\nParsing modifications...")
    modifications = parser.parse_modifications(MODIFICATIONS_CSV, "alpha_modifications.json")
    
    # Statistics
    stats = parser.get_statistics(mutations, modifications)
    print("\n" + "="*60)
    print("STATISTICS")
    print("="*60)
    print(f"Mutations: {stats['mutations']}")
    print(f"Modifications: {stats['modifications']}")
    
    # Save stats
    with open("alpha_tubulin_stats.json", "w") as f:
        json.dump(stats, f, indent=2)


if __name__ == "__main__":
    main()