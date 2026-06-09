"""Deterministic smoke test for the grounded assistant backend.

Run with Neo4j up and .env sourced:

    set -a; source .env; set +a
    python3 scripts/smoketest_assistant.py            # retrieval tools only (no LLM)
    python3 scripts/smoketest_assistant.py --llm 1JFF # also run the orchestrator loop (needs OPENAI_API_KEY)

The retrieval section is deterministic (no model) and is the fastest way to
catch Cypher / field-access mistakes in api/nl_translator/retrieval.py. Pass a
PDB id present in your DB as the optional positional arg (default 1JFF).
"""
from __future__ import annotations

import json
import sys

from api.nl_translator import retrieval as R


def show(title, obj):
    print(f"\n=== {title} ===")
    print(json.dumps(obj, indent=2, default=str)[:1500])


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    rcsb = args[0] if args else "1JFF"
    run_llm = "--llm" in sys.argv

    show("get_facets", R.run_retrieval_tool("get_facets", {}))
    show(f"get_structure_chains({rcsb})", R.run_retrieval_tool("get_structure_chains", {"rcsb_id": rcsb}))
    show("get_binding_site TA1/tubulin_beta", R.run_retrieval_tool("get_binding_site", {"chemical_id": "TA1", "family": "tubulin_beta"}))
    show("get_binding_site GTP/tubulin_alpha", R.run_retrieval_tool("get_binding_site", {"chemical_id": "GTP", "family": "tubulin_alpha"}))
    show("count_modifications tubulin_alpha (human)", R.run_retrieval_tool("count_modifications", {"family": "tubulin_alpha", "species_tax_ids": [9606]}))
    show("count_variants tubulin_beta", R.run_retrieval_tool("count_variants", {"family": "tubulin_beta"}))
    show("resolve_structure human alpha", R.run_retrieval_tool("resolve_structure", {"organism_id": 9606, "family": "tubulin_alpha"}))
    show("find_structures has TA1", R.run_retrieval_tool("find_structures", {"has_ligand_ids": ["TA1"]}))

    if run_llm:
        from api.nl_translator.orchestrator import run_assistant
        chains = R.run_retrieval_tool("get_structure_chains", {"rcsb_id": rcsb}).get("chains", [])
        ctx = {
            "page": "structure", "rcsb_id": rcsb, "view_mode": "monomer",
            "chain_ids": [c["auth_asym_id"] for c in chains],
            "active_monomer_chain": chains[0]["auth_asym_id"] if chains else None,
            "active_family": chains[0]["family"] if chains else None,
            "ligand_keys": sorted({l for c in chains for l in c["bound_ligand_chem_ids"]}),
        }
        for q in [
            "how many PTMs are recorded for human alpha-tubulin?",
            "focus the taxol binding site",
            "show me human structures with taxol",
        ]:
            res = run_assistant(q, ctx)
            print(f"\n### LLM: {q!r}")
            print(json.dumps(res.model_dump(), indent=2, default=str)[:2000])


if __name__ == "__main__":
    main()
