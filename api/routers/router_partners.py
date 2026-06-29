# api/routers/router_partners.py
from fastapi import APIRouter, HTTPException
from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import CanonicalPartnerSite

router_partners = APIRouter()


@router_partners.get(
    "/canonical-site/{map_family}/{tubulin_family}",
    response_model=CanonicalPartnerSite,
    operation_id="get_partner_canonical_site",
)
def get_partner_canonical_site(map_family: str, tubulin_family: str):
    """
    Aggregate the tubulin interface contacted by a MAP family across all structures.

    Returns per-master-alignment-position frequencies for the residues of
    `tubulin_family` that the `map_family` MAP contacts, across every structure in
    the database. Use it to paint a "canonical" MAP footprint onto any tubulin
    structure, even one where that MAP is not bound.

    Master positions are scoped to `tubulin_family` (alpha and beta have separate
    alignments). 404 if there is no data or too few structures to be canonical.
    """
    result = db_reader.get_canonical_partner_site(map_family, tubulin_family)
    if result is None:
        raise HTTPException(
            404,
            f"No canonical interface data for {map_family} on {tubulin_family}. "
            f"Either the MAP has never been observed contacting a {tubulin_family} "
            f"polymer, or too few structures contain the interface to aggregate.",
        )
    return result
