#!/usr/bin/env python3
"""
Test the new GraphQL-based profile generation
"""

import asyncio
import json
from structural_analyzer import SpatialGridGenerator

async def test_graphql_profile():
    """Test the GraphQL profile generation"""
    
    print("🔬 Testing GraphQL Profile Generation")
    print("=" * 50)
    
    analyzer = SpatialGridGenerator()
    
    try:
        # Test profile generation
        print("📡 Fetching profile via GraphQL...")
        profile = await analyzer.get_profile("6dpu")
        
        print("✅ Profile retrieved successfully!")
        
        # Check the structure/66dpu/6o
        if "entry" in profile:
            entry = profile["entry"]
            print(f"📊 Entry ID: {entry.get('rcsb_id', 'Unknown')}")
            
            # Check polymer entities
            polymer_entities = entry.get("polymer_entities", [])
            print(f"🧬 Polymer entities: {len(polymer_entities)}")
            
            # Test tubulin chain identification
            print(f"\n🔍 Identifying tubulin chains...")
            tubulin_chains = analyzer.filter_tubulin_chains(profile)
            
            if tubulin_chains:
                print(f"✅ Success! Found tubulin chains:")
                alpha_chains = [cid for cid, ctype in tubulin_chains.items() if ctype == "alpha"]
                beta_chains = [cid for cid, ctype in tubulin_chains.items() if ctype == "beta"]
                
                print(f"   Alpha chains ({len(alpha_chains)}): {alpha_chains[:10]}{'...' if len(alpha_chains) > 10 else ''}")
                print(f"   Beta chains ({len(beta_chains)}): {beta_chains[:10]}{'...' if len(beta_chains) > 10 else ''}")
                
                # Check for additional metadata
                print(f"\n📋 Additional metadata available:")
                if "rcsb_entry_info" in entry:
                    resolution = entry["rcsb_entry_info"].get("resolution_combined")
                    print(f"   Resolution: {resolution}Å" if resolution else "   Resolution: Not available")
                
                if "exptl" in entry:
                    methods = [exp.get("method", "Unknown") for exp in entry["exptl"]]
                    print(f"   Methods: {', '.join(methods)}")
                
                if "rcsb_accession_info" in entry:
                    deposit_date = entry["rcsb_accession_info"].get("deposit_date")
                    print(f"   Deposit date: {deposit_date}" if deposit_date else "   Deposit date: Not available")
                
                # Check assemblies info
                assemblies = entry.get("assemblies", [])
                print(f"   Assemblies: {len(assemblies)}")
                
                return True
            else:
                print("❌ No tubulin chains found!")
                
                # Debug: show what entities we found
                print(f"\n🔍 Available entities:")
                for i, entity in enumerate(polymer_entities[:5]):
                    polymer_info = entity.get("rcsb_polymer_entity", {})
                    desc = polymer_info.get("pdbx_description", "No description")
                    print(f"   {i}: {desc}")
                
                return False
        else:
            print("❌ Invalid profile structure - no 'entry' key found")
            print(f"🔍 Available keys: {list(profile.keys())}")
            return False
            
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

async def compare_with_old_approach():
    """Quick comparison to make sure we get the same results"""
    print(f"\n🔄 Quick comparison with old results...")
    
    # Just check that we can run the full pipeline
    analyzer = SpatialGridGenerator()
    
    try:
        # This should now use GraphQL internally
        grid_data = await analyzer.generate_grid("6dpu")
        
        print(f"✅ Full pipeline works with GraphQL!")
        print(f"   Tubulin chains: {grid_data.metadata.get('num_tubulin_chains', 'Unknown')}")
        print(f"   Protofilaments: {grid_data.metadata.get('num_protofilaments', 'Unknown')}")
        print(f"   Connections: {grid_data.metadata.get('nterm_connections', 'Unknown')}")
        
        return True
        
    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        return False

if __name__ == "__main__":
    print("🧪 Testing GraphQL Profile System")
    print("=" * 70)
    
    async def run_tests():
        # Test 1: GraphQL profile generation
        test1_ok = await test_graphql_profile()
        
        # Test 2: Full pipeline integration
        test2_ok = await compare_with_old_approach()
        
        print(f"\n🎯 Test Results:")
        print(f"   GraphQL Profile: {'✅ PASS' if test1_ok else '❌ FAIL'}")
        print(f"   Full Pipeline: {'✅ PASS' if test2_ok else '❌ FAIL'}")
        
        if test1_ok and test2_ok:
            print(f"🎉 All tests passed! GraphQL integration successful.")
        else:
            print(f"⚠️ Some tests failed - check the output above.")
    
    asyncio.run(run_tests())