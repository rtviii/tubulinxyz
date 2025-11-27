import sys
import os
import subprocess
import tempfile
from dataclasses import dataclass
from typing import List

def detect_file_format(filename: str) -> str:
    """Detect if file is PDB or mmCIF based on extension and content"""
    
    # Check extension first
    ext = os.path.splitext(filename)[1].lower()
    if ext in ['.pdb']:
        return 'pdb'
    elif ext in ['.cif', '.mmcif']:
        return 'mmcif'
    
    # Check content for ambiguous extensions
    try:
        with open(filename, 'r') as f:
            first_lines = [f.readline().strip() for _ in range(5)]
        
        # mmCIF files typically start with data_ or have loop_ early on
        content = '\n'.join(first_lines)
        if content.startswith('data_') or 'loop_' in content:
            return 'mmcif'
        elif any(line.startswith(('HEADER', 'ATOM', 'HETATM', 'REMARK')) for line in first_lines):
            return 'pdb'
    except:
        pass
    
    # Default assumption
    return 'pdb'

def convert_pdb_to_mmcif(pdb_file: str, mmcif_file: str) -> bool:
    """
    Convert PDB to mmCIF using available tools
    Returns True if successful, False otherwise
    """
    
    print(f"Converting {pdb_file} to mmCIF format...")
    
    # Method 1: Try pymol
    try:
        pymol_cmd = f"""
import pymol
pymol.cmd.load('{pdb_file}')
pymol.cmd.save('{mmcif_file}')
pymol.cmd.quit()
"""
        result = subprocess.run([
            'python', '-c', pymol_cmd
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0 and os.path.exists(mmcif_file):
            print("✓ Converted using PyMOL")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        pass
    
    # Method 2: Try gemmi
    try:
        result = subprocess.run([
            'gemmi', 'convert', pdb_file, mmcif_file
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0 and os.path.exists(mmcif_file):
            print("✓ Converted using gemmi")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        pass
    
    # Method 3: Try BioPython
    try:
        biopython_cmd = f"""
from Bio.PDB import PDBParser, MMCIFIO
parser = PDBParser(QUIET=True)
structure = parser.get_structure('structure', '{pdb_file}')
io = MMCIFIO()
io.set_structure(structure)
io.save('{mmcif_file}')
"""
        result = subprocess.run([
            'python', '-c', biopython_cmd
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0 and os.path.exists(mmcif_file):
            print("✓ Converted using BioPython")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        pass
    
    # Method 4: Try simple format conversion (basic, may not work for all cases)
    try:
        basic_conversion_cmd = f"""
# Simple PDB to mmCIF header conversion - basic fallback
with open('{pdb_file}', 'r') as f_in, open('{mmcif_file}', 'w') as f_out:
    f_out.write('data_structure\\n')
    f_out.write('#\\n')
    
    # Write basic mmCIF structure
    f_out.write('loop_\\n')
    f_out.write('_atom_site.group_PDB\\n')
    f_out.write('_atom_site.id\\n') 
    f_out.write('_atom_site.type_symbol\\n')
    f_out.write('_atom_site.label_atom_id\\n')
    f_out.write('_atom_site.label_alt_id\\n')
    f_out.write('_atom_site.label_comp_id\\n')
    f_out.write('_atom_site.label_asym_id\\n')
    f_out.write('_atom_site.label_entity_id\\n')
    f_out.write('_atom_site.label_seq_id\\n')
    f_out.write('_atom_site.pdbx_PDB_ins_code\\n')
    f_out.write('_atom_site.Cartn_x\\n')
    f_out.write('_atom_site.Cartn_y\\n')
    f_out.write('_atom_site.Cartn_z\\n')
    f_out.write('_atom_site.occupancy\\n')
    f_out.write('_atom_site.B_iso_or_equiv\\n')
    f_out.write('_atom_site.auth_seq_id\\n')
    f_out.write('_atom_site.auth_comp_id\\n')
    f_out.write('_atom_site.auth_asym_id\\n')
    f_out.write('_atom_site.auth_atom_id\\n')
    f_out.write('_atom_site.pdbx_PDB_model_num\\n')
    
    atom_id = 1
    for line in f_in:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            group_pdb = line[0:6].strip()
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()
            x = line[30:38].strip()
            y = line[38:46].strip() 
            z = line[46:54].strip()
            occupancy = line[54:60].strip() or '1.00'
            b_factor = line[60:66].strip() or '0.00'
            element = line[76:78].strip() or atom_name[0]
            
            f_out.write(f"{{group_pdb}} {{atom_id}} {{element}} {{atom_name}} . {{res_name}} {{chain_id}} 1 {{res_seq}} ? {{x}} {{y}} {{z}} {{occupancy}} {{b_factor}} {{res_seq}} {{res_name}} {{chain_id}} {{atom_name}} 1\\n")
            atom_id += 1
"""
        result = subprocess.run([
            'python', '-c', basic_conversion_cmd
        ], capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0 and os.path.exists(mmcif_file):
            print("✓ Converted using basic conversion")
            return True
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.CalledProcessError):
        pass
    
    return False

@dataclass
class ComputedResidue:
    """Represents a computed residue with metadata"""
    chain_id: str
    residue_number: int
    method: str
    confidence: float = 0.8

def add_computed_metadata_to_mmcif(mmcif_file: str, output_file: str, computed_residues: List[ComputedResidue]):
    """Add computed residue metadata to existing mmCIF file
    
    Args:
        mmcif_file: Input mmCIF file path
        output_file: Output mmCIF file path  
        computed_residues: List of ComputedResidue objects
    """
    
    # Read original mmCIF
    with open(mmcif_file, 'r') as f:
        content = f.read()
    
    # Generate metadata block
    metadata = "\n# Computed residue annotations\n"
    metadata += "loop_\n"
    metadata += "_pdbx_computed_residue.auth_asym_id\n"
    metadata += "_pdbx_computed_residue.auth_seq_id\n"
    metadata += "_pdbx_computed_residue.method\n"
    metadata += "_pdbx_computed_residue.confidence\n"
    
    # Add each computed residue
    for residue in computed_residues:
        metadata += f"{residue.chain_id}  {residue.residue_number}  '{residue.method}'  {residue.confidence}\n"
    
    metadata += "\n"
    
    # Write output file
    with open(output_file, 'w') as f:
        f.write(content)
        f.write(metadata)
    
    print(f"Added {len(computed_residues)} computed residue annotations to {output_file}")

def create_tubulin_computed_residues() -> List[ComputedResidue]:
    """Create the computed residue list for your tubulin structure"""
    
    computed_residues = []
    
    # K40 loop in alpha-tubulin (residues 38-46) - completed with Modeller
    for res_num in range(38, 47):  # 38-46 inclusive
        computed_residues.append(ComputedResidue(
            chain_id='A',
            residue_number=res_num,
            method='Modeller',
            confidence=0.7
        ))
    
    # C-terminal tail in alpha-tubulin (residues 437-451) - MD simulation
    for res_num in range(437, 452):  # 437-451 inclusive
        computed_residues.append(ComputedResidue(
            chain_id='A', 
            residue_number=res_num,
            method='MD_simulation',
            confidence=0.8
        ))
    
    # C-terminal tail in beta-tubulin (residues 427-450) - MD simulation  
    for res_num in range(427, 451):  # 427-450 inclusive
        computed_residues.append(ComputedResidue(
            chain_id='B',
            residue_number=res_num, 
            method='MD_simulation',
            confidence=0.8
        ))
    
    return computed_residues

if __name__ == "__main__":
    if len(sys.argv) not in [3, 4]:
        print("Usage: python add_metadata.py <input.pdb|input.cif> <output.cif> [--keep-temp]")
        print("\nThis will add computed residue metadata for:")
        print("- Chain A: K40 loop (38-46) via Modeller")
        print("- Chain A: C-terminal tail (437-451) via MD simulation") 
        print("- Chain B: C-terminal tail (427-450) via MD simulation")
        print("\nSupports both PDB and mmCIF input. PDB files will be converted to mmCIF first.")
        print("Use --keep-temp to keep intermediate conversion files.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    keep_temp = len(sys.argv) == 4 and sys.argv[3] == '--keep-temp'
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    # Detect input format
    input_format = detect_file_format(input_file)
    print(f"Detected input format: {input_format.upper()}")
    
    # Handle conversion if needed
    if input_format == 'pdb':
        # Create temporary mmCIF file
        temp_mmcif = tempfile.NamedTemporaryFile(suffix='.cif', delete=False)
        temp_mmcif_path = temp_mmcif.name
        temp_mmcif.close()
        
        print("Input is PDB format, converting to mmCIF...")
        if not convert_pdb_to_mmcif(input_file, temp_mmcif_path):
            print("Error: Failed to convert PDB to mmCIF")
            print("Please install one of: PyMOL, gemmi, or BioPython")
            print("Or manually convert the file to mmCIF format")
            os.unlink(temp_mmcif_path)
            sys.exit(1)
        
        mmcif_input = temp_mmcif_path
    else:
        mmcif_input = input_file
    
    try:
        # Create the tubulin computed residues
        computed_residues = create_tubulin_computed_residues()
        
        # Add to mmCIF file
        add_computed_metadata_to_mmcif(mmcif_input, output_file, computed_residues)
        
        print(f"\n✓ Successfully created {output_file}")
        print(f"\nSummary:")
        print(f"- K40 loop (Modeller): {len([r for r in computed_residues if r.method == 'Modeller'])} residues")
        print(f"- MD simulation: {len([r for r in computed_residues if r.method == 'MD_simulation'])} residues")
        print(f"- Total computed residues: {len(computed_residues)}")
        
    finally:
        # Clean up temporary file
        if input_format == 'pdb':
            if keep_temp:
                print(f"Temporary mmCIF file kept at: {temp_mmcif_path}")
            else:
                try:
                    os.unlink(temp_mmcif_path)
                except:
                    pass
