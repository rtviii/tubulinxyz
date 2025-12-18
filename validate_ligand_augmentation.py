#!/usr/bin/env python3
"""
Validate that ligand neighborhood files have correct MA index augmentation.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field

sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.etl.assets import TubulinStructureAssets, GlobalOps
from lib.types import SequenceIngestionEntry


@dataclass
class ResidueComparison:
    chain: str
    auth_seq_id: int
    master_index: int
    comp_id: str
    
    @property
    def indices_differ(self) -> bool:
        return self.auth_seq_id != self.master_index


@dataclass
class ValidationResult:
    rcsb_id: str
    ligand_file: str
    total_polymer_residues: int
    augmented_residues: int
    missing_augmentation: int
    residues_where_indices_differ: List[ResidueComparison] = field(default_factory=list)
    mismatched_vs_expected: List[Tuple[str, int, int, int]] = field(default_factory=list)
    errors: List[str] = field(default_factory=list)
    
    @property
    def is_valid(self) -> bool:
        return len(self.mismatched_vs_expected) == 0 and len(self.errors) == 0


def load_ingestion_data(ingestion_path: Path) -> Dict[str, SequenceIngestionEntry]:
    if not ingestion_path.exists():
        return {}
    with open(ingestion_path) as f:
        raw = json.load(f)
    return {
        entity_id: SequenceIngestionEntry.model_validate(entry)
        for entity_id, entry in raw.items()
    }


def build_chain_to_entity_map(profile_path: Path) -> Dict[str, str]:
    if not profile_path.exists():
        return {}
    with open(profile_path) as f:
        profile = json.load(f)
    
    chain_to_entity = {}
    for pp in profile.get("polypeptides", []):
        chain_to_entity[pp["auth_asym_id"]] = pp["entity_id"]
    return chain_to_entity


def validate_ligand_file(
    ligand_path: Path,
    chain_to_entity: Dict[str, str],
    ingestion_data: Dict[str, SequenceIngestionEntry],
) -> ValidationResult:
    rcsb_id = ligand_path.stem.split("_")[0]
    result = ValidationResult(
        rcsb_id=rcsb_id,
        ligand_file=ligand_path.name,
        total_polymer_residues=0,
        augmented_residues=0,
        missing_augmentation=0,
    )
    
    try:
        with open(ligand_path) as f:
            data = json.load(f)
        
        if isinstance(data, list):
            if not data:
                result.errors.append("Empty ligand extraction (no interactions found)")
                return result
            data = data[0]
        
        # Build reverse lookup maps
        entity_auth_to_ma: Dict[str, Dict[int, int]] = {}
        for entity_id, entry in ingestion_data.items():
            entity_auth_to_ma[entity_id] = entry.build_auth_to_ma_map()
        
        seen_residues = set()
        
        def process_residue(auth_asym_id: str, auth_seq_id: int, comp_id: str, found_ma: Optional[int]):
            nonlocal result
            
            # Skip non-polymer residues
            if auth_asym_id not in chain_to_entity:
                return
            
            key = (auth_asym_id, auth_seq_id)
            if key in seen_residues:
                return
            seen_residues.add(key)
            
            result.total_polymer_residues += 1
            entity_id = chain_to_entity[auth_asym_id]
            
            if entity_id not in entity_auth_to_ma:
                return
            
            auth_to_ma = entity_auth_to_ma[entity_id]
            expected_ma = auth_to_ma.get(auth_seq_id)
            
            if found_ma is not None:
                result.augmented_residues += 1
                
                # Check if it matches expected
                if expected_ma is not None and found_ma != expected_ma:
                    result.mismatched_vs_expected.append(
                        (auth_asym_id, auth_seq_id, found_ma, expected_ma)
                    )
                
                # Record if auth_seq_id differs from master_index
                if auth_seq_id != found_ma:
                    result.residues_where_indices_differ.append(
                        ResidueComparison(auth_asym_id, auth_seq_id, found_ma, comp_id)
                    )
            else:
                if expected_ma is not None:
                    result.missing_augmentation += 1
        
        # Process neighborhood
        for res_tuple in data.get("neighborhood", []):
            if len(res_tuple) < 3:
                continue
            auth_asym_id, auth_seq_id, comp_id = res_tuple[0], res_tuple[1], res_tuple[2]
            found_ma = res_tuple[3] if len(res_tuple) >= 4 else None
            process_residue(auth_asym_id, auth_seq_id, comp_id, found_ma)
        
        # Process interactions
        for interaction in data.get("interactions", []):
            for participant in interaction.get("participants", []):
                if len(participant) < 5:
                    continue
                auth_asym_id, auth_seq_id, comp_id = participant[0], participant[1], participant[2]
                is_ligand = participant[4]
                if is_ligand:
                    continue
                found_ma = participant[5] if len(participant) >= 6 else None
                process_residue(auth_asym_id, auth_seq_id, comp_id, found_ma)
    
    except Exception as e:
        result.errors.append(str(e))
    
    return result


def validate_structure(rcsb_id: str) -> List[ValidationResult]:
    assets = TubulinStructureAssets(rcsb_id)
    results = []
    
    ingestion_path = Path(assets.paths.sequence_ingestion)
    ingestion_data = load_ingestion_data(ingestion_path)
    
    profile_path = Path(assets.paths.profile)
    chain_to_entity = build_chain_to_entity_map(profile_path)
    
    ligand_files = assets.paths.all_ligand_neighborhoods()
    
    for lf in ligand_files:
        result = validate_ligand_file(Path(lf), chain_to_entity, ingestion_data)
        results.append(result)
    
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Validate ligand MA augmentation")
    parser.add_argument("--rcsb-id", help="Validate specific structure")
    parser.add_argument("--limit", type=int, default=100, help="Max structures to check")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--show-differences", "-d", action="store_true", 
                        help="Show all residues where auth_seq_id != master_index")
    args = parser.parse_args()
    
    if args.rcsb_id:
        structures = [args.rcsb_id.upper()]
    else:
        structures = GlobalOps.list_profiles()[:args.limit]
    
    print(f"Validating {len(structures)} structures...\n")
    
    total_files = 0
    valid_files = 0
    empty_files = 0
    files_with_mismatches = 0
    files_missing_augmentation = 0
    
    all_differences: List[Tuple[str, str, ResidueComparison]] = []
    all_mismatches: List[Tuple[str, str, Tuple]] = []
    
    for rcsb_id in structures:
        results = validate_structure(rcsb_id)
        
        for r in results:
            total_files += 1
            
            if "Empty" in str(r.errors):
                empty_files += 1
                continue
            
            if r.is_valid and r.missing_augmentation == 0:
                valid_files += 1
            
            if r.mismatched_vs_expected:
                files_with_mismatches += 1
                for m in r.mismatched_vs_expected:
                    all_mismatches.append((r.rcsb_id, r.ligand_file, m))
            
            if r.missing_augmentation > 0:
                files_missing_augmentation += 1
            
            # Collect differences
            for diff in r.residues_where_indices_differ:
                all_differences.append((r.rcsb_id, r.ligand_file, diff))
            
            if args.verbose and (r.mismatched_vs_expected or r.errors or r.residues_where_indices_differ):
                print(f"{r.rcsb_id}/{r.ligand_file}:")
                print(f"  Polymer residues: {r.total_polymer_residues}, Augmented: {r.augmented_residues}")
                
                if r.residues_where_indices_differ:
                    print(f"  Residues where auth_seq_id != MA index: {len(r.residues_where_indices_differ)}")
                
                if r.mismatched_vs_expected:
                    print(f"  VALIDATION ERRORS (found != expected): {len(r.mismatched_vs_expected)}")
                
                if r.errors:
                    print(f"  Errors: {r.errors}")
                print()
    
    # Summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Structures checked: {len(structures)}")
    print(f"Ligand files checked: {total_files}")
    print(f"Empty files (no extraction): {empty_files}")
    print(f"Valid files: {valid_files}")
    print(f"Files with MA calculation mismatches: {files_with_mismatches}")
    print(f"Files missing augmentation: {files_missing_augmentation}")
    
    if all_mismatches:
        print(f"\n** VALIDATION ERRORS (MA index != expected from ingestion): {len(all_mismatches)}")
        for rcsb_id, ligand_file, (chain, auth, found, expected) in all_mismatches[:20]:
            print(f"  {rcsb_id}/{ligand_file}: chain={chain} auth_seq_id={auth} found_MA={found} expected_MA={expected}")
    
    # Show differences where auth_seq_id != master_index
    print(f"\n" + "=" * 70)
    print(f"RESIDUES WHERE auth_seq_id != master_index: {len(all_differences)}")
    print("=" * 70)
    
    if all_differences:
        if args.show_differences or len(all_differences) <= 50:
            for rcsb_id, ligand_file, diff in all_differences:
                print(f"  {rcsb_id}/{ligand_file}: chain={diff.chain} {diff.comp_id} "
                      f"auth_seq_id={diff.auth_seq_id} -> MA_index={diff.master_index} "
                      f"(delta={diff.master_index - diff.auth_seq_id})")
        else:
            print(f"  (showing first 30, use -d to see all)")
            for rcsb_id, ligand_file, diff in all_differences[:30]:
                print(f"  {rcsb_id}/{ligand_file}: chain={diff.chain} {diff.comp_id} "
                      f"auth_seq_id={diff.auth_seq_id} -> MA_index={diff.master_index} "
                      f"(delta={diff.master_index - diff.auth_seq_id})")
        
        # Stats on deltas
        deltas = [d.master_index - d.auth_seq_id for _, _, d in all_differences]
        if deltas:
            print(f"\n  Delta stats: min={min(deltas)}, max={max(deltas)}, unique deltas={len(set(deltas))}")
    else:
        print("  None found - all auth_seq_id values equal their master_index")
        print("  (This could mean: numbering always starts at 1, or augmentation isn't working)")


if __name__ == "__main__":
    main()