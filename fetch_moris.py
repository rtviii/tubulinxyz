import requests
import time
import sys

def fetch_and_save_fasta(id_list, filename, batch_size=50):
    """
    Fetches sequences from UniProt in batches using the stream endpoint,
    appends them to the specified file, and calculates length statistics.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/stream"
    total_ids = len(id_list)
    all_lengths = []
    
    print(f"--- Processing {filename} ({total_ids} IDs) ---")
    
    # Clear/Create file initially
    with open(filename, 'w') as f:
        f.write("")
    for i in range(0, total_ids, batch_size):
        batch = id_list[i:i + batch_size]
        print(f"  Fetching batch {i//batch_size + 1} ({len(batch)} IDs)...")
        query = " OR ".join([f"accession:{pid}" for pid in batch])
        params = {
            "format": "fasta",
            "query": query
        }
        
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            
            if response.text.strip():
                # Parse FASTA text to get stats
                entries = response.text.strip().split('>')
                for entry in entries:
                    if not entry.strip(): continue
                    lines = entry.splitlines()
                    if len(lines) > 1:
                        seq = "".join(lines[1:])
                        all_lengths.append(len(seq))

                with open(filename, 'a') as f:
                    f.write(response.text)
            else:
                print("    Warning: No data returned for this batch (IDs might be obsolete).")
                
        except requests.exceptions.RequestException as e:
            print(f"    Error fetching batch: {e}")
        
        # Be polite to the API
        time.sleep(0.5)
        
    if all_lengths:
        print(f"  Stats for {filename}:")
        print(f"    - Count: {len(all_lengths)}")
        print(f"    - Min Length: {min(all_lengths)}")
        print(f"    - Max Length: {max(all_lengths)}")
    else:
        print("  No valid sequences found.")

    print(f"  Done. Saved to {filename}\n")

# ==========================================
# 1. HARDCODED ID LISTS
# ==========================================

alpha_ids = [
    "A0A5B9TE06", "A0A5B9T7A6", "A0A5B9T794", "A0A5B9T5X2", "A0A5B9T627", 
    "A0A5B9TEC1", "Q65C79", "Q65C78", "A0A4U6STT5", "A0A4U6W1W4", "A0A4U6T6L3", 
    "D5MRD7", "D5MRD8", "D5MRD9", "D5MRE0", "Q8WRT6", "P41351", "I7LY16", 
    "Q23WP5", "Q22CD2", "P06603", "P06604", "P06605", "P06606", "Q9VED6", 
    "O22347", "O22348", "O22349", "O18688", "P34690", "P91910", "G5EDD4", 
    "Q19490", "P91873", "O18154", "P52274", "Q20221", "Q71U36", "P68363", 
    "Q9BQE3", "P0DPH7", "P0DPH8", "Q6PEY2", "P68366", "Q9NY65", "A6NHL2", 
    "O22661", "O22660", "A0A059T4F2", "A0A059T4U6", "A0A059T4Y0", "A0A059T4F4", 
    "A0A059T464", "P11139", "B9DGT7", "Q56WH1", "Q0WV25", "B9DHQ0", "P29511", 
    "P28752", "Q53M52", "P68362", "P68361", "P68365", "P09204", "P09205", 
    "Q9DFR8", "Q9DD79", "Q9DFR7", "Q9DFR9", "Q78C10", "Q9YHW2", "Q9YHW1", 
    "Q9DFT3", "Q7XYS1", "Q84KQ3", "A0A125YFV3", "P04688", "P04689", "P09733", 
    "P09734", "P68369", "P05213", "P68373", "P05214", "P68368", "Q9JJZ2", 
    "A0A125YQG0", "A0A125YTP7", "A0A125YQ52"
]

beta_ids = [
    "A0A385HDR4", "P12411", "Q56YW9", "Q9ASR0", "P24636", "P29513", "P29514", 
    "P29515", "P29516", "P29517", "O17921", "P12456", "P52275", "P41937", 
    "G5EF01", "Q18817", "Q9DD57", "P04690", "P69893", "A0A125YTB7", "Q84KQ2", 
    "Q7XYS0", "Q24560", "P61857", "P08841", "Q9VAX7", "Q9VRX3", "Q9ZPP0", 
    "Q9ZPN9", "Q9ZPN8", "Q9ZPN7", "Q9N2N6", "C0L7F0", "C0L7F2", "C0L7F1", 
    "P07437", "Q9H4B7", "Q13885", "Q9BVA1", "Q13509", "P04350", "P68371", 
    "Q9BUF5", "Q3ZCM7", "A6NNZ2", "A0A5B9T5X9", "A0A5B9T5T1", "A2AQ07", 
    "Q7TMM9", "Q9CWF2", "Q9ERD7", "Q9D6F9", "P68372", "P99024", "Q922F4", 
    "P36221", "Q78C12", "Q9DFT5", "Q9DFT6", "Q43594", "Q8H7U1", "Q40665", 
    "P45960", "P46265", "Q76FS3", "P37832", "Q76FS2", "P02557", "P05219", 
    "Q65C77", "Q65C76", "A0A4U6T2L9", "A0A4V6D9T7", "A0A4U6UB67", "A0A4U6SU51", 
    "A0A4U6V2J9", "A0A4U6WIH1", "P41352", "Q24D63", "Q24D62", "I6U4Q6", 
    "Q22ED2", "A0A125YYZ6", "A0A125YJU4", "A0A125YWG5", "P10653", "P10874",
    "A0A226BMV8", "A0A226BIW8", "E2RFJ7"
]

gamma_ids = [
    "P38557", "P38558", "P34475", "Q39582", "G3HLY0", "P23257", "P42271", 
    "A3F2R1", "P23258", "Q9NRH3", "P83887", "Q8VCK3", "O49068", "P53378", 
    "P25295", "A0A4U6", "O00849", "S8F111", "P18695", "O96769", "Q9TYH0"
]

# ==========================================
# 2. EXECUTION
# ==========================================

tubulin_dict = {
    "alpha_tubulin.fasta": alpha_ids,
    "beta_tubulin.fasta":  beta_ids,
    "gamma_tubulin.fasta": gamma_ids
}

if __name__ == "__main__":
    print("Starting UniProt Sequence Fetcher...\n")
    
    for filename, id_list in tubulin_dict.items():
        fetch_and_save_fasta(id_list, filename)
        
    print("All downloads complete.")
