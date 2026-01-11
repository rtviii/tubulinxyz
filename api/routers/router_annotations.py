# api/routers/router_annotations.py
from fastapi import APIRouter, HTTPException, Query
from typing import List, Dict, Any, Optional
from neo4j import Transaction

from neo4j_tubxz.db_lib_reader import db_reader

router_annotations = APIRouter()


# =============================================================================
# Cypher Query Functions
# =============================================================================


def get_mutations_at_position(position: int, family: str) -> List[Dict[str, Any]]:
    """Get all mutations at a specific master alignment position for a family."""
    query = """
    MATCH (e:PolypeptideEntity)-[:HAS_MUTATION]->(m:Mutation)
    WHERE m.master_index = $position AND e.family = $family
    RETURN DISTINCT {
        master_index: m.master_index,
        from_residue: m.from_residue,
        to_residue: m.to_residue,
        uniprot_id: m.uniprot_id,
        species: m.species,
        phenotype: m.phenotype,
        database_source: m.database_source,
        reference_link: m.reference_link,
        keywords: m.keywords,
        rcsb_id: e.parent_rcsb_id
    } AS mutation
    """
    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            return [dict(r["mutation"]) for r in tx.run(query, {"position": position, "family": family})]
        return session.execute_read(run)


def get_mutations_for_polymer(rcsb_id: str, auth_asym_id: str) -> List[Dict[str, Any]]:
    """Get all mutations for a specific polymer chain."""
    query = """
    MATCH (i:PolypeptideInstance {parent_rcsb_id: $rcsb_id, auth_asym_id: $auth_asym_id})
    MATCH (i)-[:INSTANCE_OF]->(e:PolypeptideEntity)-[:HAS_MUTATION]->(m:Mutation)
    RETURN {
        master_index: m.master_index,
        from_residue: m.from_residue,
        to_residue: m.to_residue,
        uniprot_id: m.uniprot_id,
        species: m.species,
        phenotype: m.phenotype,
        database_source: m.database_source,
        reference_link: m.reference_link,
        keywords: m.keywords,
        rcsb_id: e.parent_rcsb_id
    } AS mutation
    ORDER BY m.master_index
    """
    with db_reader.adapter.driver.session() as session:
        def run(tx: Transaction):
            return [dict(r["mutation"]) for r in tx.run(query, {
                "rcsb_id": rcsb_id.upper(),
                "auth_asym_id": auth_asym_id
            })]
        return session.execute_read(run)


def get_interactions_at_position(position: int, family: str) -> List[Dict[str, Any]]:
    """Get ligand interactions at a master alignment position."""
    query = """
    MATCH (ni:NonpolymerInstance)-[:HAS_INTERACTION]->(ix:Interaction)-[:INTERACTS_WITH]->(pi:PolypeptideInstance)
    MATCH (pi)-[:INSTANCE_OF]->(e:PolypeptideEntity)
    WHERE ix.master_index = $position AND e.family = $family
    MATCH (ni)-[:INSTANCE_OF]->(ne:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
    RETURN DISTINCT {
        master_index: ix.master_index,
        interaction_type: ix.type,
        residue_auth_seq_id: ix.auth_seq_id,
        residue_comp_id: ix.auth_comp_id,
        atom_id: ix.atom_id,
        ligand_id: c.chemical_id,
        ligand_name: c.chemical_name,
        structure_id: ni.parent_rcsb_id,
        chain_id: pi.auth_asym_id
    } AS interaction
    """
    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return [
                dict(r["interaction"])
                for r in tx.run(query, {"position": position, "family": family})
            ]

        return session.execute_read(run)


def get_interactions_for_polymer(
    rcsb_id: str, auth_asym_id: str
) -> List[Dict[str, Any]]:
    """Get all ligand interactions for a specific polymer chain."""
    query = """
    MATCH (ni:NonpolymerInstance)-[:HAS_INTERACTION]->(ix:Interaction)-[:INTERACTS_WITH]->(pi:PolypeptideInstance)
    WHERE pi.parent_rcsb_id = $rcsb_id AND pi.auth_asym_id = $auth_asym_id
    MATCH (ni)-[:INSTANCE_OF]->(ne:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
    RETURN {
        master_index: ix.master_index,
        interaction_type: ix.type,
        residue_auth_seq_id: ix.auth_seq_id,
        residue_comp_id: ix.auth_comp_id,
        atom_id: ix.atom_id,
        ligand_id: c.chemical_id,
        ligand_name: c.chemical_name,
        ligand_chain: ni.auth_asym_id
    } AS interaction
    ORDER BY ix.master_index
    """
    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return [
                dict(r["interaction"])
                for r in tx.run(
                    query, {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id}
                )
            ]

        return session.execute_read(run)


def get_ligand_neighborhoods_for_polymer(
    rcsb_id: str, auth_asym_id: str
) -> List[Dict[str, Any]]:
    """Get all ligands in the neighborhood of a polymer chain with their nearby residues."""
    query = """
    MATCH (ni:NonpolymerInstance)-[r:NEAR_POLYMER]->(pi:PolypeptideInstance)
    WHERE pi.parent_rcsb_id = $rcsb_id AND pi.auth_asym_id = $auth_asym_id
    MATCH (ni)-[:INSTANCE_OF]->(ne:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
    RETURN {
        ligand_id: c.chemical_id,
        ligand_name: c.chemical_name,
        ligand_chain: ni.auth_asym_id,
        ligand_auth_seq_id: ni.auth_seq_id,
        nearby_residues: r.residues,
        drugbank_id: c.drugbank_id
    } AS neighborhood
    ORDER BY c.chemical_id
    """
    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return [
                dict(r["neighborhood"])
                for r in tx.run(
                    query, {"rcsb_id": rcsb_id.upper(), "auth_asym_id": auth_asym_id}
                )
            ]

        return session.execute_read(run)


def get_ligand_neighborhoods_at_position(
    position: int, family: str
) -> List[Dict[str, Any]]:
    """
    Get ligands whose neighborhood includes a specific master alignment position.
    This requires joining through interactions since NEAR_POLYMER doesn't store master_index.
    """
    query = """
    MATCH (ni:NonpolymerInstance)-[:HAS_INTERACTION]->(ix:Interaction)-[:INTERACTS_WITH]->(pi:PolypeptideInstance)
    MATCH (pi)-[:INSTANCE_OF]->(e:PolypeptideEntity)
    WHERE ix.master_index = $position AND e.family = $family
    MATCH (ni)-[:INSTANCE_OF]->(ne:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
    OPTIONAL MATCH (ni)-[r:NEAR_POLYMER]->(pi)
    RETURN DISTINCT {
        ligand_id: c.chemical_id,
        ligand_name: c.chemical_name,
        structure_id: ni.parent_rcsb_id,
        ligand_chain: ni.auth_asym_id,
        polymer_chain: pi.auth_asym_id,
        nearby_residues: r.residues,
        drugbank_id: c.drugbank_id
    } AS neighborhood
    """
    with db_reader.adapter.driver.session() as session:

        def run(tx: Transaction):
            return [
                dict(r["neighborhood"])
                for r in tx.run(query, {"position": position, "family": family})
            ]

        return session.execute_read(run)


# =============================================================================
# Position-based Endpoints
# =============================================================================


@router_annotations.get("/mutations/{family}/{position}")
async def get_mutations_at_position_endpoint(
    family: str, position: int
) -> Dict[str, Any]:
    """
    Get all mutations at a specific master alignment position.
    """
    try:
        mutations = get_mutations_at_position(position, family)
        return {
            "position": position,
            "family": family,
            "count": len(mutations),
            "mutations": mutations,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/interactions/{family}/{position}")
async def get_interactions_at_position_endpoint(
    family: str, position: int
) -> Dict[str, Any]:
    """
    Get all ligand interactions at a specific master alignment position.
    """
    try:
        interactions = get_interactions_at_position(position, family)
        return {
            "position": position,
            "family": family,
            "count": len(interactions),
            "interactions": interactions,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/neighborhoods/{family}/{position}")
async def get_neighborhoods_at_position_endpoint(
    family: str, position: int
) -> Dict[str, Any]:
    """
    Get ligand neighborhoods that include a specific master alignment position.
    """
    try:
        neighborhoods = get_ligand_neighborhoods_at_position(position, family)
        return {
            "position": position,
            "family": family,
            "count": len(neighborhoods),
            "neighborhoods": neighborhoods,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/all/{family}/{position}")
async def get_all_annotations_at_position(family: str, position: int) -> Dict[str, Any]:
    """
    Get mutations, interactions, and ligand neighborhoods at a specific position.
    """
    try:
        mutations = get_mutations_at_position(position, family)
        interactions = get_interactions_at_position(position, family)
        neighborhoods = get_ligand_neighborhoods_at_position(position, family)

        return {
            "position": position,
            "family": family,
            "mutations": mutations,
            "interactions": interactions,
            "neighborhoods": neighborhoods,
            "total_count": len(mutations) + len(interactions) + len(neighborhoods),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/range/{family}")
async def get_annotations_in_range(
    family: str,
    start: int = Query(..., description="Start position (inclusive)"),
    end: int = Query(..., description="End position (inclusive)"),
) -> Dict[str, Any]:
    """
    Get all annotations within a position range.
    """
    try:
        mut_query = """
        MATCH (e:PolypeptideEntity)-[:HAS_MUTATION]->(m:Mutation)
        WHERE m.master_index >= $start AND m.master_index <= $end AND e.family = $family
        RETURN m.master_index AS position, collect(m {
            .from_residue, .to_residue, .uniprot_id, .species, .phenotype
        }) AS mutations
        ORDER BY position
        """

        ix_query = """
        MATCH (ni:NonpolymerInstance)-[:HAS_INTERACTION]->(ix:Interaction)-[:INTERACTS_WITH]->(pi:PolypeptideInstance)
        MATCH (pi)-[:INSTANCE_OF]->(e:PolypeptideEntity)
        WHERE ix.master_index >= $start AND ix.master_index <= $end AND e.family = $family
        MATCH (ni)-[:INSTANCE_OF]->(:NonpolymerEntity)-[:DEFINED_BY_CHEMICAL]->(c:Chemical)
        RETURN ix.master_index AS position, collect({
            type: ix.type, ligand_id: c.chemical_id, structure: ni.parent_rcsb_id
        }) AS interactions
        ORDER BY position
        """

        params = {"start": start, "end": end, "family": family}

        with db_reader.adapter.driver.session() as session:

            def run(tx):
                muts = {
                    r["position"]: r["mutations"] for r in tx.run(mut_query, params)
                }
                ixs = {
                    r["position"]: r["interactions"] for r in tx.run(ix_query, params)
                }
                return muts, ixs

            mutations_by_pos, interactions_by_pos = session.execute_read(run)

        all_positions = set(mutations_by_pos.keys()) | set(interactions_by_pos.keys())
        data = {}
        for pos in sorted(all_positions):
            data[pos] = {
                "mutations": mutations_by_pos.get(pos, []),
                "interactions": interactions_by_pos.get(pos, []),
            }

        return {
            "family": family,
            "range": {"start": start, "end": end},
            "positions_with_annotations": len(data),
            "data": data,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# =============================================================================
# Polymer-based Endpoints
# =============================================================================


@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/mutations")
async def get_polymer_mutations(rcsb_id: str, auth_asym_id: str) -> Dict[str, Any]:
    """Get all mutations for a specific polymer chain."""
    try:
        mutations = get_mutations_for_polymer(rcsb_id, auth_asym_id)
        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "count": len(mutations),
            "mutations": mutations,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/interactions")
async def get_polymer_interactions(rcsb_id: str, auth_asym_id: str) -> Dict[str, Any]:
    """Get all ligand interactions for a specific polymer chain."""
    try:
        interactions = get_interactions_for_polymer(rcsb_id, auth_asym_id)
        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "count": len(interactions),
            "interactions": interactions,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/neighborhoods")
async def get_polymer_neighborhoods(rcsb_id: str, auth_asym_id: str) -> Dict[str, Any]:
    """Get all ligand neighborhoods for a specific polymer chain."""
    try:
        neighborhoods = get_ligand_neighborhoods_for_polymer(rcsb_id, auth_asym_id)
        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "count": len(neighborhoods),
            "neighborhoods": neighborhoods,
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router_annotations.get("/polymer/{rcsb_id}/{auth_asym_id}/all")
async def get_polymer_all_annotations(
    rcsb_id: str, auth_asym_id: str
) -> Dict[str, Any]:
    """Get mutations, interactions, and ligand neighborhoods for a polymer chain."""
    try:
        mutations = get_mutations_for_polymer(rcsb_id, auth_asym_id)
        interactions = get_interactions_for_polymer(rcsb_id, auth_asym_id)
        neighborhoods = get_ligand_neighborhoods_for_polymer(rcsb_id, auth_asym_id)

        return {
            "rcsb_id": rcsb_id.upper(),
            "auth_asym_id": auth_asym_id,
            "mutations": mutations,
            "interactions": interactions,
            "neighborhoods": neighborhoods,
            "total_count": len(mutations) + len(interactions) + len(neighborhoods),
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
