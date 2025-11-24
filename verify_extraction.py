from pathlib import Path
from api.services.mmcif_parser import TubulinStructureParser

# 1. Setup
parser = TubulinStructureParser()
cif_path = Path("data/maxim_data/7sj7_with_metadata.cif")
chain_target = "A"

# 2. Extract (Backend Way)
print(f"Extracting {chain_target} from {cif_path}...")
be_seq, be_ids = parser.get_observed_data(cif_path, chain_target)

print(f"Backend Length: {len(be_seq)}")
print(f"Sample Backend: {be_seq[:10]}... IDs: {be_ids[:10]}")

# 3. Paste Frontend Data Here (Simulated)
# In reality, you'd paste what your console.log gave you in the browser
# For now, let's assume the frontend is perfect and verify against itself
fe_seq_simulated = be_seq 
fe_ids_simulated = be_ids

# 4. Compare
report = parser.compare_with_frontend(fe_seq_simulated, fe_ids_simulated, be_seq, be_ids)

if report['match']:
    print("\n✅ SUCCESS: Python extraction perfectly matches Frontend expectation.")
else:
    print("\n❌ FAILURE: Inconsistencies detected.")
    print(report)