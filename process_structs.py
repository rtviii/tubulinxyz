from Bio.PDB.MMCIFParser import  MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBList import PDBList
import os

def extract_chain_from_structure(file_path, chain_id, output_path):
    """
    Extract a specific chain from a PDB/CIF structure file and save it.
    
    Parameters:
    file_path (str): Path to the input structure file
    chain_id (str): Chain identifier to extract (e.g., 'A')
    output_path (str): Path where the extracted chain will be saved
    """
    
    # Determine file format and use appropriate parser
    if file_path.lower().endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    
    # Parse the structure
    structure_id = os.path.splitext(os.path.basename(file_path))[0]
    structure = parser.get_structure(structure_id, file_path)
    
    # Create PDBIO object for writing
    io = PDBIO()
    
    # Find and extract the specified chain
    chain_found = False
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                print(f"Found chain {chain_id} in {structure_id}")
                print(f"Chain {chain_id} has {len(chain)} residues")
                
                # Create a new structure with just this chain
                io.set_structure(chain)
                io.save(output_path)
                print(f"Chain {chain_id} saved to {output_path}")
                chain_found = True
                break
        if chain_found:
            break
    
    if not chain_found:
        print(f"Chain {chain_id} not found in {structure_id}")
        # List available chains
        available_chains = []
        for model in structure:
            for chain in model:
                available_chains.append(chain.id)
        print(f"Available chains: {available_chains}")

def main():
    """
    Extract chain A from both tubulin structures
    """
    # Structure files (assuming they're in current directory)
    structures = [
        ("6S8M.cif", "6S8M_chainA.pdb"),
        ("5JCO.cif", "5JCO_chainA.pdb")
    ]
    
    chain_id = "A"
    
    for input_file, output_file in structures:
        if os.path.exists(input_file):
            print(f"\nProcessing {input_file}...")
            extract_chain_from_structure(input_file, chain_id, output_file)
        else:
            print(f"\nWarning: {input_file} not found in current directory")
            print("Please make sure the CIF files are downloaded and in the same directory as this script")

if __name__ == "__main__":
    main()

# Alternative: If you want to download the structures directly from PDB
def download_and_extract_chains():
    """
    Download structures from PDB and extract chain A
    """
    
    # Download structures
    pdbl = PDBList()
    pdb_codes = ["6S8M", "5JCO"]
    
    for pdb_code in pdb_codes:
        print(f"\nDownloading {pdb_code}...")
        # Download as mmCIF format
        pdbl.retrieve_pdb_file(pdb_code, file_format='mmCif', pdir='.')
        
        # The downloaded file will have a specific naming convention
        cif_file = f"{pdb_code.lower()}.cif"
        output_file = f"{pdb_code}_chainA.pdb"
        
        if os.path.exists(cif_file):
            extract_chain_from_structure(cif_file, "A", output_file)
