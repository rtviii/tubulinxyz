from typing_extensions import Literal
import typing
from Bio.SeqRecord import SeqRecord
from pydantic import BaseModel
import os
from ete3 import NCBITaxa

import threading

from etl.constants import NCBI_TAXA_SQLITE



_thread_local = threading.local()

def get_ncbi():
    """Get a thread-local NCBITaxa instance"""
    if not hasattr(_thread_local, 'ncbi'):
        _thread_local.ncbi = NCBITaxa(dbfile=NCBI_TAXA_SQLITE)
    return _thread_local.ncbi

TAXID_BACTERIA = 2
TAXID_EUKARYOTA = 2759
TAXID_ARCHAEA = 2157
PhylogenyRank = Literal[
    "superkingdom",
    "realm",
    'cellular root',
    'acellular root',
    'serogroup',
    "phylum",
    "domain",
    "class",
    "order",
    "clade",
    "family",
    "genus",
    "species",
    "strain",
    "isolate",
    "subspecies",
    "no rank",
    "suborder",
    "kingdom",
    "subfamily",
    "subgenus",
    "subphylum",
    "infraorder",
    "superorder",
    "superclass",
    "superfamily",
    "parvorder",
    "cohort",
    "infraclass",
    "subclass",
    "subkingdom",
    "species group",
    "tribe",
    "species subgroup",
    "subcohort",
    "subtribe",
]


# def find_n_closest_relatives(
#     _taxids_base: list[int], taxid_target: str, n_neighbors: int = 10
# ) -> list[int]:
#     """Given a set of taxids and a target taxid, return a list of the [n_neighbors] phylogenetically closest to the target."""

#     taxids_base = list(
#         map(lambda x: str(x), _taxids_base)
#     )  # Ensure all taxids are strings because that's what ete's ncbi expects

#     tree = ncbi.get_topology(list(set([*taxids_base, str(taxid_target)])))
#     target_node = tree.search_nodes(name=str(taxid_target))[0]
#     phylo_all_nodes = [
#         (other_node.name, tree.get_distance(target_node, other_node))
#         for other_node in tree.traverse()
#     ]

#     phylo_extant_nodes = filter(lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)
#     phylo_sorted_nodes = sorted(
#         phylo_extant_nodes, key=lambda x: x[1]
#     )  # Sort by phylogenetic distance (tuples are (taxid, phylo_dist) ex. ('4932', 3.0))
#     nbr_taxids = list(map(lambda tax_phydist: tax_phydist[0], phylo_sorted_nodes))

#     if len(nbr_taxids) < n_neighbors:
#         return list(map(lambda x: int(x), nbr_taxids[1:]))
#     else:
#         return list(map(lambda x: int(x), nbr_taxids[1 : n_neighbors + 1]))


def find_n_closest_relatives(_taxids_base: list[int], taxid_target: str, n_neighbors: int = 10) -> list[int]:
    """
    Find the N phylogenetically closest organisms to a target organism from a given set of candidates.
    
    This function constructs a phylogenetic tree containing both the target taxid and all base taxids,
    then finds the N closest organisms to the target based on tree distance. The target itself is 
    excluded from the results.
    
    Args:
        _taxids_base (list[int]): List of candidate taxids to search through. These are the potential
            neighbors from which the closest ones will be selected.
        taxid_target (str): The reference taxid to compare against. This is the organism whose
            neighborhood we want to explore.
        n_neighbors (int, optional): Number of neighbors to return. Defaults to 10.
            If len(_taxids_base) < n_neighbors, returns all available neighbors.
    
    Returns:
        list[int]: List of the N phylogenetically closest taxids to the target, sorted by increasing
            phylogenetic distance. The length of the list will be min(n_neighbors, len(_taxids_base)).
    
    Example:
        >>> base_taxids = [9598, 9601, 10090, 10116]  # Chimp, Orangutan, Mouse, Rat
        >>> target = "9606"  # Human
        >>> find_n_closest_relatives(base_taxids, target, 2)
        [9598, 9601]  # Returns Chimp and Orangutan as closest to Human
    
    Notes:
        - Uses NCBI taxonomy database via ete3 for phylogenetic relationships
        - Distance is based on the topology of the phylogenetic tree
        - The target taxid is always excluded from the results, even if it exists in _taxids_base
    """
    # Convert all taxids to strings for ete3 compatibility
    taxids_base = list(map(lambda x: str(x), _taxids_base))


    ncbi = get_ncbi()
    # Get topology tree containing target and all base taxids
    tree = ncbi.get_topology(list(set([*taxids_base, str(taxid_target)])))
    target_node = tree.search_nodes(name=str(taxid_target))[0]
    
    # Calculate distances from target to all other nodes
    phylo_all_nodes = [(other_node.name, tree.get_distance(target_node, other_node)) 
                       for other_node in tree.traverse()]

    # Filter to keep only nodes that were in our original base set
    phylo_extant_nodes = filter(lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)
    
    # Sort by phylogenetic distance
    phylo_sorted_nodes = sorted(phylo_extant_nodes, key=lambda x: x[1])  # Sort by phylogenetic distance
    nbr_taxids = list(map(lambda tax_phydist: tax_phydist[0], phylo_sorted_nodes))

    # Return either all available neighbors or n_neighbors, whichever is smaller
    if len(nbr_taxids) < n_neighbors:
        return list(map(lambda x: int(x), nbr_taxids[1:]))
    else:
        return list(map(lambda x: int(x), nbr_taxids[1:n_neighbors + 1]))


def get_lineage_distance(source_taxid: int, target_taxid: int) -> float:
    """
    Calculate distance between two taxa based on their complete taxonomic lineages.
    
    The distance is computed by:
    1. Getting full lineages for both taxa
    2. Finding their most recent common ancestor
    3. Counting the number of distinct ranks between each taxon and their common ancestor
    
    Args:
        source_taxid (int): First taxid
        target_taxid (int): Second taxid
    
    Returns:
        float: Distance based on taxonomic lineage divergence. Higher numbers indicate
              taxa that diverged earlier in evolutionary history.
    
    Example:
        >>> # Compare E. coli (562) to Human (9606)
        >>> dist = get_lineage_distance(562, 9606)
        >>> print(f"Distance: {dist}")
        >>> # Will show a much larger distance than before since they diverged very early
    """
    # Get complete lineages for both taxa
    source_lineage = set(Taxid.get_lineage(source_taxid))
    target_lineage = set(Taxid.get_lineage(target_taxid))
    
    # Find shared lineage (all common ancestors)
    shared_lineage = source_lineage.intersection(target_lineage)
    
    # Get unique ranks in each lineage (steps since divergence)
    source_unique = source_lineage - shared_lineage
    target_unique = target_lineage - shared_lineage
    
    # The total evolutionary distance is the sum of unique steps in both lineages
    # This represents how many taxonomic ranks have diverged since their last common ancestor
    distance = len(source_unique) + len(target_unique)
    
    return float(distance)

def print_lineage_comparison(source_taxid: int, target_taxid: int):
    """
    Helper function to visualize the lineage comparison between two taxa.
    
    Args:
        source_taxid (int): First taxid
        target_taxid (int): Second taxid
    """
    source_lineage = Taxid.get_lineage(source_taxid)
    target_lineage = Taxid.get_lineage(target_taxid)
    
    # Convert to sets for intersection
    source_set = set(source_lineage)
    target_set = set(target_lineage)
    shared = source_set.intersection(target_set)
    
    print(f"\nComparing {Taxid.get_name(source_taxid)} to {Taxid.get_name(target_taxid)}:")
    
    print("\nShared lineage (common ancestors):")
    for taxid in sorted(shared, key=lambda x: source_lineage.index(x)):
        print(f"  {Taxid.get_name(taxid)} ({Taxid.rank(taxid)})")
    
    print(f"\nUnique to {Taxid.get_name(source_taxid)}:")
    for taxid in sorted(source_set - shared, key=lambda x: source_lineage.index(x)):
        print(f"  {Taxid.get_name(taxid)} ({Taxid.rank(taxid)})")
    
    print(f"\nUnique to {Taxid.get_name(target_taxid)}:")
    for taxid in sorted(target_set - shared, key=lambda x: target_lineage.index(x)):
        print(f"  {Taxid.get_name(taxid)} ({Taxid.rank(taxid)})")
    
    distance = get_lineage_distance(source_taxid, target_taxid)
    print(f"\nLineage distance: {distance}")



class Taxid:
    @staticmethod
    def is_descendant_of(parent_taxid: int, target_taxid: int) -> bool:
        ncbi = get_ncbi()
        lineage = ncbi.get_lineage(target_taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        return False if parent_taxid not in lineage else True


    @staticmethod
    def get_phylogenetic_distance(source_taxid: int, target_taxid: int) -> float:
        """
        Calculate the phylogenetic distance between two taxids in the NCBI taxonomy tree.
        
        Args:
            source_taxid (int): First taxid
            target_taxid (int): Second taxid
        
        Returns:
            float: Phylogenetic distance between the two taxids. Distance represents the number of 
                  taxonomic ranks that need to be traversed to get from one taxid to another through 
                  their most recent common ancestor.
        
        Raises:
            LookupError: If either taxid is invalid in NCBI database
        
        Example:
            >>> dist = get_phylogenetic_distance(9606, 9598)  # Human vs Chimp
            >>> print(f"Distance: {dist:.4f}")
            Distance: 0.0132
        """
        # Convert taxids to strings for ete3
        str_source = str(source_taxid)
        str_target = str(target_taxid)
        
        ncbi = get_ncbi()
        # Get topology tree containing both taxids
        tree = ncbi.get_topology([str_source, str_target])
        
        source_node = tree.search_nodes(name=str_source)[0]
        target_node = tree.search_nodes(name=str_target)[0]
        print(source_node, target_node)
        
        # Calculate distance
        return tree.get_distance(source_node, target_node)


    @staticmethod
    def get_name(taxid):


        ncbi = get_ncbi()
        return list(ncbi.get_taxid_translator([taxid]).values())[0]

    @staticmethod
    def get_lineage(
        taxid, include_only: None | list[PhylogenyRank] = None
    ) -> list[int]:
        """Return ncbi lineage, except filter out the ranks that are not among the @PhylogenyRank."""
        # lin = list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )
        # lin = list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )
        ncbi = get_ncbi()
        lin = ncbi.get_lineage(taxid)
        if include_only is not None:
            return list(filter(lambda x: Taxid.rank(x) in include_only, lin))
        return lin if lin is not None else []

    @staticmethod
    def rank(taxid: int) -> PhylogenyRank:
        """Given a @taxid, return the rank of the taxid"""

        ncbi = get_ncbi()
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage)[taxid]

    @staticmethod
    def coerce_to_rank(taxid: int, target_rank: PhylogenyRank) -> int | None:
        """Given a @taxid and a @rank, return the taxid of the first ancestor of @taxid that is at @rank"""


        ncbi = get_ncbi()
        lineage = ncbi.get_lineage(taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        for item in lineage:
            rank = ncbi.get_rank([item])[item]
            if rank == target_rank:
                return item

        raise IndexError("Taxid {} does not have a {} level".format(taxid, target_rank))

    @staticmethod
    def superkingdom(
        taxid: int,
    ) -> typing.Literal["bacteria", "eukaryota", "archaea", "virus"]:
        match (
            Taxid.is_descendant_of(TAXID_EUKARYOTA, taxid),
            Taxid.is_descendant_of(TAXID_BACTERIA, taxid),
            Taxid.is_descendant_of(TAXID_ARCHAEA, taxid),
        ):
            case (False, False, True):
                return "archaea"
            case (False, True, False):
                return "bacteria"
            case (True, False, False):
                return "eukaryota"
            case (False, False, False):
                print("Probably a virus")
                return "virus"
            case _:
                raise ValueError(
                    "Taxid {} is not a descendant of any of the three domains".format(
                        taxid
                    )
                )

    @staticmethod
    def coerce_all_to_rank(taxids: list[int], level: PhylogenyRank) -> list[int]:
        """Given a list of taxids, return a list of the same taxids but coerced to the given lineage level(rank)."""
        new = []
        for taxid in taxids:
            try:
                new.append(Taxid.coerce_to_rank(taxid, level))
            except Exception as e:
                print(e)
                raise Exception(e)
        return new

    @staticmethod
    def get_descendants_of(parent: int, targets: list[int]):
        """Given a @parent taxid and a list of @taxids, return the subset of @taxids that are descendants of @parent"""
        descendants = set()
        ncbi = get_ncbi()
        for tax_id in targets:

            lineage = ncbi.get_lineage(tax_id)
            if lineage == None:
                raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
            if parent in lineage:
                descendants.add(tax_id)
        return descendants


class PhylogenyNode(BaseModel):

    @staticmethod
    def from_taxid(taxid: int):
        return PhylogenyNode(
            ncbi_tax_id=taxid,
            scientific_name=Taxid.get_name(taxid),
            rank=Taxid.rank(taxid),
        )

    def __hash__(self) -> int:
        return self.ncbi_tax_id

    def get_lineage(self) -> list[int]:
        return Taxid.get_lineage(self.ncbi_tax_id)

    def get_lineage_nodes(self) -> list[int]:
        _ = []
        for taxid in self.get_lineage():
            if Taxid.rank(taxid) not in typing.get_args(PhylogenyRank):
                continue
            _.append(self.from_taxid(taxid))
        return _

    ncbi_tax_id: int
    scientific_name: str
    rank: PhylogenyRank


# ? ----------------------------------------------{ Subcomponent Types }------------------------------------------------
