#!/usr/bin/env python3
"""
Quick test with improved parameters
"""

import asyncio
from structural_analyzer import SpatialGridGenerator

async def quick_test():
    analyzer = SpatialGridGenerator()
    
    print("🔧 Testing with improved parameters:")
    print("   - N-terminus radius: 7Å (was 5Å)")  
    print("   - Min contacts: 2 (was 3)")
    print("   - Better start chain detection")
    print("   - Smarter isolated chain grouping")
    print()
    
    try:
        grid_data = await analyzer.generate_grid("6o2t")
        
        print(f"\n📊 IMPROVED RESULTS:")
        print(f"   Connections: {grid_data.metadata['nterm_connections']}/104 ({grid_data.metadata['nterm_connections']/104*100:.1f}%)")
        print(f"   Protofilaments: {grid_data.metadata['num_protofilaments']} (target: ~13)")
        
        # Count protofilament sizes
        pf_sizes = {}
        for subunit in grid_data.subunits:
            pf_idx = subunit.protofilament
            pf_sizes[pf_idx] = pf_sizes.get(pf_idx, 0) + 1
        
        size_counts = {}
        for size in pf_sizes.values():
            size_counts[size] = size_counts.get(size, 0) + 1
        
        print(f"\n🧵 Protofilament size distribution:")
        for size, count in sorted(size_counts.items(), reverse=True):
            print(f"   Size {size}: {count} protofilaments")
        
        # Success metrics
        long_pfs = sum(1 for size in pf_sizes.values() if size >= 4)
        connection_rate = grid_data.metadata['nterm_connections'] / 104
        
        print(f"\n🎯 Success metrics:")
        print(f"   Protofilaments ≥4 chains: {long_pfs}")
        print(f"   Connection rate: {connection_rate:.1%}")
        
        if long_pfs >= 10 and connection_rate > 0.5:
            print(f"✅ Much better! Ready for PyMOL validation")
        elif connection_rate > 0.4:
            print(f"👍 Improved, may need further tuning")
        else:
            print(f"⚠️ Still needs work - consider 8Å radius or contact threshold = 1")
            
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    asyncio.run(quick_test())