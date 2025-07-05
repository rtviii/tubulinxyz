#!/usr/bin/env python3
"""
Test the fixed grid layout with layer-based clustering
"""

import asyncio
from structural_analyzer import SpatialGridGenerator

# RCSB_ID = "6O2T"  # Example PDB ID for testing
# RCSB_ID = "4TV9"  # Example PDB ID for testing
RCSB_ID = "6DPU"  # Example PDB ID for testing
RCSB_ID = "6wvr"  # Example PDB ID for testing
RCSB_ID = "8qv0"  # Example PDB ID for testing
RCSB_ID = "6fkj"  # Example PDB ID for testing
RCSB_ID = "6u0h"  # Example PDB ID for testing
RCSB_ID = "8vrk"  # Example PDB ID for testing

async def test_fixed_grid():
    """Test the layer-based Z-clustering fix"""
    print("🔧 Testing FIXED Grid Layout (Layer-based clustering)")
    print("=" * 70)
    analyzer = SpatialGridGenerator()
    
    try:
        grid_data = await analyzer.generate_grid(RCSB_ID)
        
        print(f"✅ Grid generation successful!")
        
        # Analyze the improved grid layout
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
        
        print(f"\n📊 IMPROVED Grid Layout:")
        print(f"   Grid dimensions: {max_pf + 1} × {max_su + 1}")
        print(f"   Total positions: {len(grid_positions)}")
        print(f"   Expected positions: {(max_pf + 1) * (max_su + 1)}")
        
        print(f"\n🏗️ Layer distribution:")
        for layer in sorted(layer_counts.keys()):
            count = layer_counts[layer]
            print(f"   Layer {layer}: {count} subunits")
        
        # Check if we have a proper 2D grid now
        if max_su > 0:
            print(f"✅ SUCCESS! Multi-layer grid created ({max_su + 1} layers)")
            
            # Show sample positions from different layers
            print(f"\n📋 Sample grid positions by layer:")
            for layer in range(min(3, max_su + 1)):  # Show first 3 layers
                layer_positions = [(pf, su) for (pf, su) in grid_positions.keys() if su == layer]
                if layer_positions:
                    sample_pos = sorted(layer_positions)[:5]  # First 5 protofilaments
                    print(f"   Layer {layer}: positions {sample_pos}")
                    
                    # Show what's at the first position in this layer
                    if sample_pos:
                        first_pos = sample_pos[0]
                        subunits_here = grid_positions[first_pos]
                        chain_ids = [s.auth_asym_id for s in subunits_here]
                        types = [s.monomerType[0] for s in subunits_here]
                        print(f"     Position {first_pos}: {'/'.join(types)} {chain_ids}")
            
            # Calculate grid quality
            expected_total = (max_pf + 1) * (max_su + 1)
            actual_total = len(grid_positions)
            quality = (actual_total / expected_total) * 100
            
            print(f"\n🎯 Grid Quality:")
            print(f"   Occupancy: {actual_total}/{expected_total} ({quality:.1f}%)")
            
            if quality > 90:
                print(f"   ✅ Excellent! Nearly complete grid")
            elif quality > 75:
                print(f"   👍 Good grid with some gaps")
            else:
                print(f"   ⚠️ Sparse grid, may need tuning")
                
        else:
            print(f"❌ Still flat! All subunits at layer 0")
            print(f"   This suggests Z-coordinates aren't being properly differentiated")
        
        return grid_data
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    asyncio.run(test_fixed_grid())