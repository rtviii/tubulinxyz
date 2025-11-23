import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from api.services.structure_parser import TubulinStructureParser

parser = TubulinStructureParser()
pdb = "1JFF"
chain = "A"

print(f"Fetching {pdb} from RCSB...")
path = parser.fetch_cif_to_temp(pdb)

if path:
    print(f"Downloaded to {path}")
    seq, ids = parser.get_observed_data(path, chain)
    print(f"Parsed Sequence Length: {len(seq)}")
    print(f"First 5 residues: {seq[:5]}")
    print(f"First 5 IDs: {ids[:5]}")
    
    import os
    os.unlink(path)
    print("Temp file cleaned up.")
else:
    print("Download failed.")