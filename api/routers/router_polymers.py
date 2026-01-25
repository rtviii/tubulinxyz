# api/routers/router_polymers.py
from fastapi import APIRouter, Query, HTTPException
from typing import Optional, List

from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import PolypeptideEntityFilters, PolypeptideListResponse

router_polymers = APIRouter()


@router_polymers.get("", response_model=PolypeptideListResponse, operation_id="list_polymers")
def list_polymers(
    cursor: Optional[str] = Query(None),
    limit: int = Query(25, ge=1, le=500),
    # Structure-level
    resolution_min: Optional[float] = Query(None, alias="resMin"),
    resolution_max: Optional[float] = Query(None, alias="resMax"),
    year_min: Optional[int] = Query(None, alias="yearMin"),
    year_max: Optional[int] = Query(None, alias="yearMax"),
    source_organism_ids: Optional[List[int]] = Query(None, alias="sourceTaxa"),
    # Entity-level
    family: Optional[List[str]] = Query(None),
    uniprot_accession: Optional[str] = Query(None, alias="uniprot"),
    sequence_contains: Optional[str] = Query(None, alias="motif"),
    sequence_length_min: Optional[int] = Query(None, alias="seqLenMin"),
    sequence_length_max: Optional[int] = Query(None, alias="seqLenMax"),
    has_variants: Optional[bool] = Query(None, alias="hasVariants"),
):
    """List polypeptide entities with filters and pagination."""
    filters = PolypeptideEntityFilters(
        cursor=cursor,
        limit=limit,
        resolution_min=resolution_min,
        resolution_max=resolution_max,
        year_min=year_min,
        year_max=year_max,
        source_organism_ids=source_organism_ids,
        family=family,
        uniprot_accession=uniprot_accession,
        sequence_contains=sequence_contains,
        sequence_length_min=sequence_length_min,
        sequence_length_max=sequence_length_max,
        has_variants=has_variants,
    )
    return db_reader.list_polypeptide_entities(filters)
