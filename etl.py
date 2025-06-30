

import asyncio
import requests

from gql_querystrings import AssemblyIdentificationString, EntryInfoString, NonpolymerEntitiesString, PolymerEntitiesString, WholeStructureTemplate


gql_string = """
{
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
}"""

import json
import os
import requests

class TubulinETLCollector:
    """
    Simple ETL Pipeline for tubulin structures - just dumps raw GraphQL data to JSON
    """
    
    rcsb_id: str
    
    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id
    
    def query_rcsb_api(self, gql_string: str) -> dict:
        """Query the RCSB GraphQL API"""
        reqstring = "https://data.rcsb.org/graphql?query={}".format(gql_string)
        _resp = requests.get(reqstring)
        resp = _resp.json()
        
        if "data" in resp:
            return resp["data"]
        else:
            raise Exception("No data found for query: {}".format(gql_string))
    
    async def generate_profile(self, output_dir: str = "./profiles") -> dict:
        """Generate raw structure profile using WholeStructureTemplate"""
        
        # Query using the WholeStructureTemplate
        raw_data = self.query_rcsb_api(
            WholeStructureTemplate.replace("$RCSB_ID", self.rcsb_id)
        )
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Save raw JSON
        output_path = os.path.join(output_dir, f"{self.rcsb_id}_profile.json")
        with open(output_path, "w") as outfile:
            json.dump(raw_data, outfile, indent=4)
        
        return raw_data

async def generate_profile(rcsb_id: str) -> dict:
    """Generate profile for tubulin structure"""
    profile = await TubulinETLCollector(rcsb_id).generate_profile()
    return profile

print(asyncio.run(generate_profile("6O2T")))  # Example RCSB ID for testing
