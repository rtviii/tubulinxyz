"""Deterministic smoke test for the grounded assistant backend.

Run with Neo4j up and .env sourced:

    set -a; source .env; set +a
    python3 scripts/smoketest_assistant.py            # strip selftest + retrieval tools (no LLM)
    python3 scripts/smoketest_assistant.py --strip-only  # JUST the strip selftest (no running DB)
    python3 scripts/smoketest_assistant.py --llm 1JFF # also run the orchestrator loop (needs OPENAI_API_KEY)

The retrieval section is deterministic (no model) and is the fastest way to
catch Cypher / field-access mistakes in api/nl_translator/retrieval.py. Pass a
PDB id present in your DB as the optional positional arg (default 1JFF). The
strip selftest (`--strip-only`) needs no running Neo4j, but still needs .env
sourced — importing the orchestrator instantiates the (lazy, unconnected) driver.
"""
from __future__ import annotations

import json
import sys

from api.nl_translator import retrieval as R


def show(title, obj):
    print(f"\n=== {title} ===")
    print(json.dumps(obj, indent=2, default=str)[:1500])


def strip_leak_selftest() -> bool:
    """Deterministic (no DB, no LLM) check of orchestrator._strip_tool_call_leak —
    the sanitizer now applied to EVERY respond terminal's answer_markdown/summary,
    not just the bare-text fallback. Mirrors the two field-test symptoms: a
    tool-call fragment leaked into the prose (A1) and literal \\n escapes + a stray
    trailing quote (A2). Run alone with `--strip-only` (needs .env sourced for the
    import chain, but no running Neo4j)."""
    from api.nl_translator.orchestrator import _strip_tool_call_leak as strip

    BN = "\\n"  # the two-character sequence the model wrongly emits for a newline
    cases = [
        ("A1 <parameter> card leak",
         'The GTP binding site spans residues 140-180.", '
         '<parameter name="cards">[{"action":"open_expert","rcsb_id":"6S8L"}]',
         lambda o: o == "The GTP binding site spans residues 140-180."),
        ("A1 viewer_actions leak",
         'Colchicine binds at the intra-dimer interface.", '
         '"viewer_actions": [{"type":"FocusBindingSite","args":{"chemical_id":"LOC"}}]',
         lambda o: o == "Colchicine binds at the intra-dimer interface."),
        ("A2 literal \\n + stray trailing quote",
         "The cryo-EM structures are:" + BN + BN + "- 9WDA" + BN + '- 9WDB"',
         lambda o: "\n" in o and o.endswith("9WDB")),
        ("double-encoded JSON string literal",
         '"Human alpha-tubulin has 482 PTM records across 101 positions.'
         + BN + BN + 'Most are phosphorylation."',
         lambda o: not o.startswith('"') and not o.endswith('"') and "482" in o),
        ("clean text passes through unchanged",
         "Beta tubulin has 12 disease variants near the taxol site.",
         lambda o: o == "Beta tubulin has 12 disease variants near the taxol site."),
    ]
    print("\n=== strip_leak_selftest (deterministic, no DB/LLM) ===")
    all_ok = True
    for label, raw, ok_fn in cases:
        out = strip(raw)
        # Universal invariant on top of the per-case predicate: no leaked marker
        # and no literal \n/\t escape may survive.
        leaked = any(m in out for m in ("<parameter", "<invoke", "<function", "<antml")) \
            or BN in out or "\\t" in out
        ok = bool(ok_fn(out)) and not leaked
        print(f"  [{'ok' if ok else 'FAIL'}] {label}")
        if not ok:
            print(f"        in : {raw!r}")
            print(f"        out: {out!r}")
            all_ok = False
    return all_ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    rcsb = args[0] if args else "1JFF"
    run_llm = "--llm" in sys.argv

    if not strip_leak_selftest():
        print("\nstrip_leak_selftest FAILED — fix _strip_tool_call_leak before trusting respond output.")
        sys.exit(1)
    if "--strip-only" in sys.argv:
        return

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
