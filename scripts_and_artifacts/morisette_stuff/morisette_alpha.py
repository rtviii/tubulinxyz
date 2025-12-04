import requests
import time

# Alpha-tubulin accession IDs from the paper
alpha_tubulins = [
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
    "P09734", "P68369", "P05213", "P68373", "P05214", "P68368", "O9IIZ2",
    "A0A125YOG0", "A0A125YTP7", "A0A125YO52"
]

def fetch_uniprot_sequence(accession):
    """Fetch a sequence from UniProt in FASTA format"""
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        print(f"Failed to fetch {accession}: {response.status_code}")
        return None

# def main():
#     output_file = "alpha_tubulins.fasta"
    
#     with open(output_file, "w") as f:
#         for i, accession in enumerate(alpha_tubulins, 1):
#             print(f"Fetching {i}/{len(alpha_tubulins)}: {accession}")
#             sequence = fetch_uniprot_sequence(accession)
            
#             if sequence:
#                 f.write(sequence)
            
#             # Be nice to UniProt's servers
#             time.sleep(0.5)
    
#     print(f"\nDone! Sequences saved to {output_file}")

def parse_fasta(filename):
    """Parse FASTA file and return list of (header, sequence) tuples"""
    sequences = []
    current_header = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header:
            sequences.append((current_header, ''.join(current_seq)))
    
    return sequences

def main():
    sequences = parse_fasta("alpha_tubulins.fasta")
    
    print(f"Total sequences: {len(sequences)}\n")
    
    # Find longest sequence
    longest_seq = max(sequences, key=lambda x: len(x[1]))
    accession = longest_seq[0].split('|')[1] if '|' in longest_seq[0] else longest_seq[0].split()[0]
    
    print(f"Longest sequence: {accession}")
    print(f"Length: {len(longest_seq[1])} amino acids")
    print(f"Full header: {longest_seq[0]}")
    print()
    
    # Show all sorted by length
    print("All sequences sorted by length:")
    print("-" * 60)
    sorted_seqs = sorted(sequences, key=lambda x: len(x[1]), reverse=True)
    for i, (header, seq) in enumerate(sorted_seqs, 1):
        acc = header.split('|')[1] if '|' in header else header.split()[0]
        print(f"{i}. {acc}: {len(seq)} aa")

if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()