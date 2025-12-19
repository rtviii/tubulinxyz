import json
import os
import typing
import requests
from typing import Any

from lib.etl.constants import TUBETL_DATA
from lib.etl.libtax import PhylogenyNode, Taxid
from lib.types import TubulinStructure 

import os
import json
from lib.etl.constants import TUBETL_DATA

class TubulinStructureAssetPaths:
    """Manages file paths for a given tubulin structure asset."""
    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()

    @property
    def base_dir(self) -> str:
        """The main directory for this structure's assets."""
        return os.path.join(TUBETL_DATA, self.rcsb_id)

    @property
    def profile(self) -> str:
        """Path to the main Pydantic JSON profile."""
        return os.path.join(self.base_dir, f"{self.rcsb_id}.json")

    @property
    def cif(self) -> str:
        """Path to the original CIF file."""
        return os.path.join(self.base_dir, f"{self.rcsb_id}.cif")

    @property
    def chains_dir(self) -> str:
        """Directory for split chain files."""
        return os.path.join(self.base_dir, "CHAINS")
    
    def ligand_neighborhood(self, comp_id: str, auth_asym_id: str) -> str:
        """Path to a specific ligand neighborhood report."""
        return os.path.join(self.base_dir, f"{self.rcsb_id}_{comp_id}_{auth_asym_id}.json")

    def all_ligand_neighborhoods(self) -> list[str]:
            """List only ligand neighborhood report files, excluding metadata."""
            import glob
            # Specifically look for the [ID]_[COMP]_[ASYM].json pattern
            pattern = os.path.join(self.base_dir, f"{self.rcsb_id}_*_*.json")
            
            excluded_suffixes = [
                f"{self.rcsb_id}.json", 
                "_classification_report.json", 
                "sequence_ingestion.json"
            ]
            
            return [
                p for p in glob.glob(pattern) 
                if not any(p.endswith(s) for s in excluded_suffixes)
            ]

    @property
    def classification_report(self) -> str:
        """Path to the HMM classification report."""
        return os.path.join(self.base_dir, f"{self.rcsb_id}_classification_report.json")

    @property
    def sequence_ingestion(self) -> str:
        """Path to the sequence ingestion results JSON."""
        return os.path.join(self.base_dir, "sequence_ingestion.json")

class TubulinStructureAssets:
    """Manager for accessing and verifying tubulin structure assets."""
    rcsb_id: str
    paths: TubulinStructureAssetPaths

    def __init__(self, rcsb_id: str) -> None:
        if not TUBETL_DATA:
            raise EnvironmentError("TUBETL_DATA environment variable is not set.")
        self.rcsb_id = rcsb_id.upper()
        self.paths = TubulinStructureAssetPaths(self.rcsb_id)

    def profile(self) -> TubulinStructure:
        """Loads the Pydantic model from the JSON profile."""
        if not os.path.exists(self.paths.profile):
            raise FileNotFoundError(f"Profile not found for {self.rcsb_id} at {self.paths.profile}")
        
        with open(self.paths.profile, "r") as f:
            return TubulinStructure.model_validate(json.load(f))

    def _verify_dir_exists(self):
        """Ensures the base directory for this structure exists."""
        if not os.path.exists(self.paths.base_dir):
            os.umask(0)
            os.makedirs(self.paths.base_dir, 0o755, exist_ok=True)
            print(f"Created asset directory at: {self.paths.base_dir}")
# --------------------------------------------------------------------------
# GQL DATA-RETRIEVAL STRINGS (from gql_querystrings.py)
# --------------------------------------------------------------------------

EntryInfoString = """
{
  entry(entry_id: "$RCSB_ID") {
  rcsb_id
  rcsb_accession_info{
    deposit_date
  }
    struct_keywords {
      pdbx_keywords
      text
    }
    rcsb_entry_info {
      resolution_combined
    }
    rcsb_external_references {
      link
      type
      id
    }
    exptl {
      method
    }
    citation {
      rcsb_authors
      year
      title
      pdbx_database_id_DOI
    }
    struct_keywords {
      pdbx_keywords
      text
    }
  }
}
"""

AssemblyIdentificationString = """
{
  entry(entry_id: "$RCSB_ID") {
    assemblies {
      rcsb_id
      nonpolymer_entity_instances {
        rcsb_nonpolymer_entity_instance_container_identifiers {
          comp_id
          auth_asym_id
          rcsb_id
          auth_seq_id
          entity_id
        }
      }
      polymer_entity_instances {
        rcsb_polymer_entity_instance_container_identifiers {
          asym_id
          auth_asym_id
          entry_id
          entity_id
        }
      }
    }
  }
}
"""

PolymerEntitiesString = """{
  entry(entry_id: "$RCSB_ID") {
    
    rcsb_id
    assemblies {
      polymer_entity_instances {
        rcsb_id
      }
    }
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        asym_ids
        auth_asym_ids
        entry_id
        entity_id
      }
      pfams {
        rcsb_pfam_accession
        rcsb_pfam_comment
        rcsb_pfam_description
      }
      rcsb_entity_source_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      rcsb_entity_host_organism {
        ncbi_taxonomy_id
        scientific_name
      }
      uniprots {
        rcsb_id
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      entity_poly {
        pdbx_seq_one_letter_code
        pdbx_seq_one_letter_code_can
        pdbx_strand_id
        rcsb_entity_polymer_type
        rcsb_sample_sequence_length
        type
      }
      rcsb_polymer_entity_annotation {
        annotation_id
        assignment_version
        description
        name
        provenance_source
        type
      }
    }
  }
}
"""

NonpolymerEntitiesString = """{
  entry(entry_id: "$RCSB_ID") {
    nonpolymer_entities {
      pdbx_entity_nonpoly {
        comp_id
        name
        entity_id
      }
      rcsb_nonpolymer_entity {
        details
        formula_weight
        pdbx_description
        pdbx_number_of_molecules
      }
      rcsb_nonpolymer_entity_container_identifiers {
        entity_id
        asym_ids
        auth_asym_ids
      }
      nonpolymer_comp {
        chem_comp {
          id
          name
          three_letter_code
        }
        rcsb_chem_comp_target {
          interaction_type
          name
          provenance_source
          reference_database_accession_code
          reference_database_name
        }
        drugbank {
          drugbank_container_identifiers {
            drugbank_id
          }
          drugbank_info {
            cas_number
            description
            indication
            mechanism_of_action
            name
            pharmacology
          }
        }
      }
    }
  }
}
"""
LigandsChemInfo = """{
  chem_comps(comp_ids: $COMP_IDS) {
    chem_comp{
      id
    }
    rcsb_chem_comp_descriptor {
      SMILES
      InChI
      InChIKey
      SMILES_stereo
    }
  }
}
"""

# --------------------------------------------------------------------------
# GQL SEARCH QUERY (from gql_querystrings.py, as a dict)
# --------------------------------------------------------------------------

TUBULIN_SEARCH_QUERY = {
    "query": {
        "type": "group",
        "logical_operator": "or",
        "nodes": [
            {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.annotation_lineage.id",
                            "operator": "exact_match",
                            "negation": False,
                            "value": "IPR036525"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.type",
                            "operator": "exact_match",
                            "value": "InterPro",
                            "negation": False
                        }
                    }
                ],
                "label": "nested-attribute"
            },
            {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.annotation_lineage.id",
                            "operator": "exact_match",
                            "negation": False,
                            "value": "IPR002452"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.type",
                            "operator": "exact_match",
                            "value": "InterPro",
                            "negation": False
                        }
                    }
                ],
                "label": "nested-attribute"
            },
            {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.annotation_lineage.id",
                            "operator": "exact_match",
                            "negation": False,
                            "value": "IPR013838"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.type",
                            "operator": "exact_match",
                            "value": "InterPro",
                            "negation": False
                        }
                    }
                ],
                "label": "nested-attribute"
            },
            {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.annotation_lineage.id",
                            "operator": "exact_match",
                            "negation": False,
                            "value": "IPR023123"
                        }
                    },
                    {
                        "type": "terminal",
                        "service": "text",
                        "parameters": {
                            "attribute": "rcsb_uniprot_annotation.type",
                            "operator": "exact_match",
                            "value": "InterPro",
                            "negation": False
                        }
                    }
                ],
                "label": "nested-attribute"
            }
        ],
        "label": "text"
    },
    "return_type": "entry",
    "request_options": {
        "return_all_hits": True,
        "results_verbosity": "compact",
        "results_content_type": ["experimental"]
    }
}

# --------------------------------------------------------------------------
# GLOBAL OPS CLASS (from global_ops.py)
# --------------------------------------------------------------------------

class GlobalOps:
    """Global utilities for managing the tubulin structure dataset."""

    @staticmethod
    def missing_profiles() -> list[str]:
        """Return a list of structures that are in the RCSB but not in the local database."""
        try:
            rcsb_structs = set(GlobalOps.current_rcsb_structs())
        except Exception as e:
            print(f"Error fetching current RCSB structs: {e}")
            return []
            
        local_profiles = set(GlobalOps.list_profiles())
        return list(rcsb_structs - local_profiles)

    @staticmethod
    def missing_db_entries(db_entries: list[str]) -> list[str]:
        """Return a list of structures that are in the RCSB but not in the local database."""
        try:
            rcsb_structs = set(GlobalOps.current_rcsb_structs())
        except Exception as e:
            print(f"Error fetching current RCSB structs: {e}")
            return []
            
        return list(rcsb_structs - set(db_entries))

    @staticmethod
    def current_rcsb_structs() -> list[str]:
        """
        Return all tubulin structures in the RCSB based on the
        project-defined search query.
        """
        rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query"
        
        # TUBE-FIX: The error was here. 
        # We pass the dict directly, not json.loads(dict)
        query_payload = TUBULIN_SEARCH_QUERY

        try:
            response = requests.post(rcsb_search_api, json=query_payload)
            response.raise_for_status() # Raise an error for bad status codes
            
            resp_json = response.json()
            
            if "result_set" not in resp_json:
                raise Exception(f"Unexpected API response: {resp_json}")

            # TUBE-UPDATE: Parse the result_set format correctly
            return sorted(resp_json["result_set"])
        
        except requests.exceptions.RequestException as e:
            print(f"Failed to query RCSB Search API: {e}")
            raise e
        except json.JSONDecodeError:
            print(f"Failed to decode JSON response from RCSB Search API: {response.text}")
            raise

    @staticmethod
    def list_profiles() -> list[str]:
        """
        List all structures that have a valid profile.json in the TUBETL_DATA dir.
        """

        if not os.path.exists(TUBETL_DATA):
            print(f"Warning: TUBETL_DATA directory not found at: {TUBETL_DATA}")
            return []
            
        profiles_exist = []
        for rcsb_id in os.listdir(TUBETL_DATA):
            if len(rcsb_id) != 4 or not os.path.isdir(os.path.join(TUBETL_DATA, rcsb_id)):
                continue
                
            profile_path = TubulinStructureAssets(rcsb_id).paths.profile
            if os.path.exists(profile_path):
                profiles_exist.append(rcsb_id)
                
        return profiles_exist

    @staticmethod
    def collect_all_taxa() -> set[PhylogenyNode]:
        """
        Scans all local profiles and collects all unique organism taxIDs
        to seed the phylogeny tree.
        """
        _ = set()
        print(f"Collecting taxa from {len(GlobalOps.list_profiles())} profiles...")
        
        for struct in GlobalOps.list_profiles():
            # TUBE-UPDATE: Use TubulinStructureAssets
            try:
                profile = TubulinStructureAssets(struct).profile()
            except Exception as e:
                print(f"Warning: Could not load profile {struct}. Skipping. Error: {e}")
                continue
            
            all_orgs = profile.src_organism_ids + profile.host_organism_ids
            
            for org_id in all_orgs:
                if org_id is None:
                    continue
                try:
                    pn = PhylogenyNode(
                        ncbi_tax_id=org_id,
                        scientific_name=Taxid.get_name(org_id),
                        rank=Taxid.rank(org_id),
                    )
                    _.add(pn)
                except Exception as e:
                    print(f"Error processing taxid {org_id} for {struct}: {e}")
                    # This can happen for obsolete or unclassified taxIDs
        
        print(f"Found {len(_)} unique taxa.")
        return _

# --------------------------------------------------------------------------
# RCSB DATA API QUERY HELPER (from queries.py)
# --------------------------------------------------------------------------

RCSB_GQL_URL = "https://data.rcsb.org/graphql"

def query_rcsb_api(gql_string: str) -> dict[str, Any]:
    """
    Submits a GraphQL query to the RCSB DATA API.
    This is different from the SEARCH API used in GlobalOps.
    
    Raises:
        Exception: If the query fails or returns no data.
    """
    try:
        response = requests.get(f"{RCSB_GQL_URL}?query={gql_string}")
        response.raise_for_status()  # Raise an exception for bad status codes
        
        resp_json = response.json()
        
        if "data" in resp_json and resp_json["data"]:
            return resp_json["data"]
        elif "errors" in resp_json:
            raise Exception(f"RCSB API returned errors: {resp_json['errors']}")
        else:
            raise Exception("No data found for query.")
            
    except requests.exceptions.RequestException as e:
        raise Exception(f"Failed to query RCSB API: {e}")
    except json.JSONDecodeError:
        raise Exception(f"Failed to decode JSON response: {response.text}")