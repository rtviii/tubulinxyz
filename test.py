#!/usr/bin/env python3
"""
Simple test script for the improved N-terminus connectivity algorithm
"""

import asyncio
import json
from pathlib import Path
from structural_analyzer import SpatialGridGenerator


async def test_connectivity(pdb_id: str = "6o2t"):
    """Test the improved connectivity algorithm"""
    
    print(f"🧪 Testing N-terminus connectivity for {pdb_id.upper()}")
    print("=" * 60)
    
    # Create analyzer and run grid generation
    analyzer = SpatialGridGenerator()
    
    try:
        grid_data = await analyzer.generate_grid(pdb_id)
        
        print(f"\n✅ SUCCESS! Generated grid for {pdb_id.upper()}")
        print(f"📊 Results:")
        print(f"   Total subunits: {len(grid_data.subunits)}")
        print(f"   Tubulin chains: {grid_data.metadata['num_tubulin_chains']}")
        print(f"   N-term connections: {grid_data.metadata['nterm_connections']}")
        print(f"   Protofilaments: {grid_data.metadata['num_protofilaments']}")
        print(f"   Structure type: {grid_data.structure_type}")
        
        # Analyze protofilament distribution
        pf_counts = {}
        alpha_counts = {}
        beta_counts = {}
        
        for subunit in grid_data.subunits:
            pf_idx = subunit.protofilament
            if pf_idx not in pf_counts:
                pf_counts[pf_idx] = 0
                alpha_counts[pf_idx] = 0
                beta_counts[pf_idx] = 0
            
            pf_counts[pf_idx] += 1
            if subunit.monomerType == "alpha":
                alpha_counts[pf_idx] += 1
            else:
                beta_counts[pf_idx] += 1
        
        print(f"\n🧵 Protofilament breakdown:")
        for pf_idx in sorted(pf_counts.keys()):
            total = pf_counts[pf_idx]
            alpha = alpha_counts[pf_idx]
            beta = beta_counts[pf_idx]
            print(f"   PF{pf_idx}: {total} subunits ({alpha}α, {beta}β)")
        
        # Check for expected microtubule structure
        expected_pfs = 13  # Typical microtubule has 13 protofilaments
        if len(pf_counts) == expected_pfs:
            print(f"🎯 Perfect! Found expected {expected_pfs} protofilaments")
        elif len(pf_counts) < expected_pfs:
            print(f"⚠️  Found {len(pf_counts)} protofilaments, expected ~{expected_pfs}")
            print(f"   Some protofilaments may not be fully connected")
        else:
            print(f"🤔 Found {len(pf_counts)} protofilaments, more than expected {expected_pfs}")
        
        # Connection efficiency
        connection_rate = grid_data.metadata['nterm_connections'] / grid_data.metadata['num_tubulin_chains']
        print(f"\n📈 Connection efficiency: {connection_rate:.1%}")
        
        if connection_rate > 0.8:
            print(f"✅ Excellent connection rate!")
        elif connection_rate > 0.6:
            print(f"👍 Good connection rate")
        else:
            print(f"⚠️  Low connection rate - may need algorithm tuning")
        
        # Check debug files
        debug_dir = Path("debug_output")
        debug_files = list(debug_dir.glob("*"))
        
        print(f"\n🔍 Debug files created: {len(debug_files)}")
        for file in debug_files:
            print(f"   - {file.name}")
        
        # PyMOL instructions
        print(f"\n🎨 To visualize in PyMOL:")
        print(f"   cd debug_output/")
        print(f"   wget https://files.rcsb.org/download/{pdb_id.upper()}.cif")
        print(f"   ./view_{pdb_id}.sh")
        print(f"   (or manually: pymol {pdb_id.upper()}.cif -d \"@{pdb_id}_visualize.pml\")")
        
        return grid_data
        
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


async def compare_old_vs_new(pdb_id: str = "6o2t"):
    """Compare results and show improvements"""
    print(f"🔬 Analysis summary for {pdb_id.upper()}")
    print("=" * 60)
    
    result = await test_connectivity(pdb_id)
    
    if result:
        print(f"\n💡 Key improvements in this version:")
        print(f"   ✅ Follows N-terminus connections to build complete protofilaments")
        print(f"   ✅ No more 'lost' monomers in the linking process")
        print(f"   ✅ Proper graph traversal instead of buggy linked list")
        print(f"   ✅ PyMOL visualization for biological validation")
        print(f"   ✅ Debug output to understand what's happening")
        
        print(f"\n🎯 Next steps:")
        print(f"   1. Check PyMOL visualization to validate 5Å N-terminus cutoff")
        print(f"   2. Verify that connected chains are biologically reasonable")
        print(f"   3. If results look good, proceed with Z-axis alignment")
        print(f"   4. If not, adjust cutoff distance or minimum contact threshold")


if __name__ == "__main__":
    import sys
    
    # Get PDB ID from command line or use default
    pdb_id = sys.argv[1] if len(sys.argv) > 1 else "6o2t"
    
    print(f"🚀 Running improved connectivity test...")
    
    # Ensure debug directory exists
    Path("debug_output").mkdir(exist_ok=True)
    
    # Run the test
    asyncio.run(compare_old_vs_new(pdb_id))