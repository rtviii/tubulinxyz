#!/usr/bin/env python3
"""
CLI tool for tubulin structure analysis
"""

import asyncio
import argparse
import sys
from pathlib import Path

# Add project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from tubulin_analyzer import SpatialGridGenerator


async def analyze_structure(pdb_id: str, debug: bool = False):
    """Analyze a PDB structure and generate grid"""
    print(f"Analyzing structure {pdb_id.upper()}")
    print("=" * 50)
    
    analyzer = SpatialGridGenerator()
    
    try:
        grid_data = await analyzer.generate_grid(pdb_id)
        
        print(f"Analysis successful!")
        
        # Analyze the grid layout
        grid_positions = {}
        layer_counts = {}
        
        for subunit in grid_data.subunits:
            pf_idx = subunit.protofilament
            su_idx = subunit.subunitIndex
            
            # Track grid positions
            pos = (pf_idx, su_idx)
            if pos not in grid_positions:
                grid_positions[pos] = []
            grid_positions[pos].append(subunit)
            
            # Count layers
            if su_idx not in layer_counts:
                layer_counts[su_idx] = 0
            layer_counts[su_idx] += 1
        
        max_pf = max(s.protofilament for s in grid_data.subunits) if grid_data.subunits else 0
        max_su = max(s.subunitIndex for s in grid_data.subunits) if grid_data.subunits else 0
        
        print(f"\nGrid Layout Analysis:")
        print(f"   Grid dimensions: {max_pf + 1} × {max_su + 1}")
        print(f"   Total positions: {len(grid_positions)}")
        print(f"   Expected positions: {(max_pf + 1) * (max_su + 1)}")
        print(f"   Structure type: {grid_data.structure_type}")
        
        if debug:
            print(f"\nLayer distribution:")
            for layer in sorted(layer_counts.keys()):
                count = layer_counts[layer]
                print(f"   Layer {layer}: {count} subunits")
            
            # Show sample positions from different layers
            print(f"\nSample grid positions by layer:")
            for layer in range(min(3, max_su + 1)):
                layer_positions = [(pf, su) for (pf, su) in grid_positions.keys() if su == layer]
                if layer_positions:
                    sample_pos = sorted(layer_positions)[:5]
                    print(f"   Layer {layer}: positions {sample_pos}")
        
        # Calculate grid quality
        expected_total = (max_pf + 1) * (max_su + 1)
        actual_total = len(grid_positions)
        quality = (actual_total / expected_total) * 100 if expected_total > 0 else 0
        
        print(f"\nGrid Quality:")
        print(f"   Occupancy: {actual_total}/{expected_total} ({quality:.1f}%)")
        
        if quality > 90:
            print(f"   Status: Excellent! Nearly complete grid")
        elif quality > 75:
            print(f"   Status: Good grid with some gaps")
        else:
            print(f"   Status: Sparse grid, may need tuning")
        
        # Show debug files created
        debug_dir = Path("debug_output")
        if debug_dir.exists():
            debug_files = list(debug_dir.glob(f"*{pdb_id}*"))
            if debug_files:
                print(f"\nDebug files created in debug_output/:")
                for file in debug_files:
                    print(f"   {file.name}")
        
        return grid_data
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        return None


def main():
    parser = argparse.ArgumentParser(description="Analyze tubulin structures and generate 2D grids")
    parser.add_argument("pdb_id", help="PDB ID to analyze (e.g., 6o2t)")
    parser.add_argument("--debug", "-d", action="store_true", help="Enable debug output")
    
    args = parser.parse_args()
    
    # Ensure debug output directory exists
    Path("debug_output").mkdir(exist_ok=True)
    
    # Run the analysis
    result = asyncio.run(analyze_structure(args.pdb_id, args.debug))
    
    if result is None:
        sys.exit(1)
    else:
        print(f"\nAnalysis completed successfully!")


if __name__ == "__main__":
    main()