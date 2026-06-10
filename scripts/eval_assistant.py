"""Behavioral eval harness for the grounded assistant.

This is the regression net the assistant has been missing. `smoketest_assistant.py`
checks the deterministic retrieval tools (Cypher / field access); THIS harness
exercises the full model loop (`run_assistant`) against a fixed case set and makes
*structural* assertions on the `AssistantResult` — never exact-string matches, since
the model is non-deterministic. What it asserts:

  - kind            the terminal (respond / clarify / cannot)
  - tools           retrieval tools that MUST appear in trace[] (grounding proof:
                    it looked the fact up instead of reciting it). A nested list is
                    an any-of group.
  - number          a digit appears in the answer/data (count questions)
  - table           data.table is populated (columns + rows)
  - action_types    viewer/suggested action types present (nested list = any-of)
  - cards_action    at least one routing card with the given action
  - query_has       a catalogue query carries ALL the named structure-filter keys

`soft` expectations warn but never fail — for model-discretion niceties (a suggested
chip, a table) we'd like but won't gate on.

One invariant is checked on every case regardless of `expect`: `no_tool_leak` — the
answer text must never carry a leaked tool-call fragment or literal \n/\t escape.

Run with Neo4j up and .env sourced:

    set -a; source .env; set +a
    python3 scripts/eval_assistant.py              # full run, 1x each
    python3 scripts/eval_assistant.py --repeat 3   # 3x each -> flake rate
    python3 scripts/eval_assistant.py 1JFF         # use a specific structure for ctx

Exit code is non-zero if any hard assertion failed (CI-friendly). Every (case, run)
is logged in full to notes/eval_runs/<ts>.jsonl for post-mortem.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
import traceback
from typing import Any, Dict, List, Optional

# Make `api` / `neo4j_tubxz` importable no matter how this script is launched
# (`python scripts/eval_assistant.py` sets sys.path[0] to scripts/, not the root).
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from api.nl_translator import retrieval as R


# ---------------------------------------------------------------------------
# Page-context presets — what the frontend would send for each page.
# ---------------------------------------------------------------------------

def _structure_ctx(rcsb: str, prefer_family: Optional[str] = None) -> Dict[str, Any]:
    """Build a structure/monomer page_context from a real structure in the DB.
    AddAnnotationTrack / AlignChain require monomer view, so we use it."""
    chains = R.run_retrieval_tool("get_structure_chains", {"rcsb_id": rcsb}).get("chains", [])
    if not chains:
        raise RuntimeError(f"{rcsb} has no chains in the DB — pick a structure that exists.")
    active = chains[0]
    if prefer_family:
        active = next((c for c in chains if c.get("family") == prefer_family), chains[0])
    return {
        "page": "structure",
        "rcsb_id": rcsb,
        "view_mode": "monomer",
        "chain_ids": [c["auth_asym_id"] for c in chains],
        "active_monomer_chain": active["auth_asym_id"],
        "active_family": active.get("family"),
        "ligand_keys": sorted({l for c in chains for l in c.get("bound_ligand_chem_ids", [])}),
        "loaded_tracks": [],
    }


def build_context(case: Dict[str, Any], rcsb: str) -> Dict[str, Any]:
    page = case["page"]
    if page == "landing":
        return {"page": "landing"}
    if page == "structure_alpha":
        return _structure_ctx(rcsb, prefer_family="tubulin_alpha")
    return _structure_ctx(rcsb)


# ---------------------------------------------------------------------------
# Cases — the 7 triage questions (notes/assistant_roadmap.md) + 3 smoketest
# carryovers. Each hits a different website<->db<->llm surface.
# ---------------------------------------------------------------------------

CASES: List[Dict[str, Any]] = [
    dict(id="ptms_human_alpha", page="structure",
         query="What PTMs are in human alpha tubulin?",
         expect=dict(kind="respond", tools=["count_modifications"], number=True),
         soft=dict(table=True, action_types=[["AddAnnotationTrack"]])),

    dict(id="colchicine_binds", page="structure",
         query="Show me where colchicine binds",
         expect=dict(kind="respond",
                     tools=[["get_binding_site", "get_binding_contacts",
                             "resolve_structure", "get_structure_chains"]]),
         soft=dict(action_types=[["AddAnnotationTrack", "FocusBindingSite"]])),

    dict(id="cryoem_human_taxol", page="landing",
         query="Find cryo-EM human structures with taxol better than 3 Å",
         expect=dict(kind="respond", tools=["find_structures"],
                     cards_action=["open_catalogue"]),
         soft=dict(query_has=["has_ligand_ids", "source_organism_ids",
                              "exp_method", "resolution_max"])),

    dict(id="variants_near_taxol", page="structure",
         query="Are there disease variants in beta tubulin near the taxol site?",
         expect=dict(kind="respond", tools=["count_variants"]),
         soft=dict(tools=[["get_binding_site", "get_binding_contacts"]])),

    dict(id="compare_ptms_human_toxo", page="structure",
         query="Compare PTMs between human and Toxoplasma alpha tubulin",
         expect=dict(kind="respond", tools=["count_modifications"]),
         soft=dict(table=True)),

    dict(id="add_toxo_to_alignment", page="structure_alpha",
         query="Add a Toxoplasma alpha-tubulin sequence to the alignment",
         expect=dict(kind="respond", tools=["resolve_structure"],
                     action_types=["AlignChain"]),
         soft=dict()),

    dict(id="kd_taxol_honest", page="structure",
         query="What is the binding affinity (Kd) of taxol?",
         expect=dict(kind="cannot"),
         soft=dict()),

    # --- smoketest carryovers ---
    dict(id="smoke_ptm_count", page="structure",
         query="how many PTMs are recorded for human alpha-tubulin?",
         expect=dict(kind="respond", tools=["count_modifications"], number=True),
         soft=dict()),

    dict(id="smoke_focus_taxol", page="structure",
         query="focus the taxol binding site",
         expect=dict(kind="respond"),
         soft=dict(action_types=[["FocusBindingSite", "AddAnnotationTrack"]])),

    dict(id="smoke_human_taxol_list", page="landing",
         query="show me human structures with taxol",
         expect=dict(kind="respond", tools=["find_structures"]),
         soft=dict(cards_action=["open_catalogue", "open_structure"])),
]


# ---------------------------------------------------------------------------
# Assertion evaluation — returns a list of failure messages (empty == pass).
# ---------------------------------------------------------------------------

def _trace_tools(result) -> List[str]:
    return [t.tool for t in result.trace]


def _action_types(result) -> List[str]:
    types = [a.type for a in result.viewer_actions]
    types += [s.action.type for s in result.suggested_actions]
    return types


def _result_text(result) -> str:
    return " ".join([
        result.answer_markdown or "",
        result.summary or "",
        json.dumps(result.data or {}, default=str),
    ])


def _require_membership(label: str, reqs, present: List[str]) -> List[str]:
    """Each req is a string (must be present) or a list (any-of must be present)."""
    fails = []
    for req in reqs:
        if isinstance(req, (list, tuple)):
            if not any(x in present for x in req):
                fails.append(f"{label}: expected any of {list(req)} in {present}")
        elif req not in present:
            fails.append(f"{label}: expected {req!r} in {present}")
    return fails


def evaluate(result, expect: Dict[str, Any]) -> List[str]:
    fails: List[str] = []

    if "kind" in expect and result.kind != expect["kind"]:
        why = result.reason or result.clarification or ""
        fails.append(f"kind: expected {expect['kind']!r}, got {result.kind!r} ({why})")

    fails += _require_membership("tools", expect.get("tools", []), _trace_tools(result))
    fails += _require_membership("action_types", expect.get("action_types", []), _action_types(result))

    if expect.get("number") and not re.search(r"\d", _result_text(result)):
        fails.append("number: expected a digit in the answer/data")

    if expect.get("table"):
        tbl = (result.data or {}).get("table") if result.data else None
        if not (tbl and tbl.get("columns") and tbl.get("rows")):
            fails.append("table: expected data.table with columns + rows")

    for ca in expect.get("cards_action", []):
        if not any(c.action == ca for c in result.cards):
            fails.append(f"cards_action: expected a card with action={ca!r}")

    qkeys = expect.get("query_has")
    if qkeys:
        ok = False
        for q in result.queries:
            fs = getattr(q, "filters_structures", None)
            if fs is None:
                continue
            d = fs.model_dump(exclude_none=True)
            if all(d.get(k) not in (None, [], "") for k in qkeys):
                ok = True
                break
        if not ok:
            fails.append(f"query_has: no catalogue query carries all of {qkeys}")

    return fails


def no_tool_leak(result) -> List[str]:
    """Hard invariant on EVERY case (case-independent, so it lives outside the
    per-case `expect`): a respond terminal's answer_markdown/summary must never
    carry a tool-call leak (`<parameter`/`<invoke`/...) or literal `\\n`/`\\t`
    escapes. Guards the sanitizer in orchestrator._build_terminal_result."""
    fails: List[str] = []
    for text, label in [(result.answer_markdown, "answer_markdown"),
                        (result.summary, "summary")]:
        if not text:
            continue
        if any(m in text for m in ("<parameter", "<invoke", "<function", "<antml")):
            fails.append(f"{label}: tool-call marker fragment leaked into text")
        if "\\n" in text or "\\t" in text:
            fails.append(f"{label}: literal \\n/\\t escape in text")
    return fails


# ---------------------------------------------------------------------------
# Deterministic retrieval pre-check — fail fast on a DB/Cypher problem before
# spending any LLM calls.
# ---------------------------------------------------------------------------

def retrieval_precheck(rcsb: str) -> bool:
    print("=== retrieval pre-check (deterministic, no model) ===")
    checks = [
        ("get_facets", {}, lambda r: bool(r.get("tubulin_families"))),
        ("get_structure_chains", {"rcsb_id": rcsb}, lambda r: bool(r.get("chains"))),
        ("count_modifications", {"family": "tubulin_alpha", "species_tax_ids": [9606]},
         lambda r: "total_records" in r),
        ("resolve_structure", {"organism_id": 9606, "family": "tubulin_alpha"},
         lambda r: isinstance(r, dict)),
    ]
    all_ok = True
    for name, args, ok_fn in checks:
        try:
            r = R.run_retrieval_tool(name, args)
            ok = ok_fn(r)
        except Exception as e:  # noqa: BLE001
            ok = False
            r = {"error": f"{type(e).__name__}: {e}"}
        print(f"  [{'ok' if ok else 'FAIL'}] {name}")
        if not ok:
            print(f"        -> {json.dumps(r, default=str)[:300]}")
            all_ok = False
    return all_ok


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def _is_fatal_provider_error(exc: Exception) -> bool:
    """Auth/quota errors (bad key, monthly limit) aren't assistant failures —
    detect them so the run aborts with a clear message instead of marking every
    remaining case as a regression."""
    s = f"{type(exc).__name__}: {exc}".lower()
    return any(k in s for k in (
        "permissiondenied", "authenticationerror", "insufficient_quota",
        "limit exceeded", "invalid api key", "code: 401", "code: 403",
    ))


def main() -> int:
    argv = sys.argv[1:]
    repeat = 1
    positional: List[str] = []
    skip_next = False
    for j, a in enumerate(argv):
        if skip_next:
            skip_next = False
            continue
        if a == "--repeat":
            if j + 1 < len(argv):
                repeat = max(1, int(argv[j + 1]))
                skip_next = True  # consume the value, don't treat it as positional
            continue
        if a.startswith("--"):
            continue
        positional.append(a)
    rcsb = positional[0] if positional else "1JFF"

    if not retrieval_precheck(rcsb):
        print("\nretrieval pre-check failed — fix the DB/Cypher layer before the LLM eval.")
        return 2

    try:
        from api.nl_translator.orchestrator import run_assistant
    except Exception as e:  # noqa: BLE001
        print(f"\nCannot import orchestrator ({e}). Skipping LLM cases.")
        return 2

    run_dir = os.path.join("notes", "eval_runs")
    os.makedirs(run_dir, exist_ok=True)
    artifact = os.path.join(run_dir, f"{int(time.time())}.jsonl")

    print(f"\n=== LLM behavioral eval ({len(CASES)} cases x {repeat}) — ctx structure {rcsb} ===")
    print(f"    artifact: {artifact}\n")

    case_pass = 0
    case_flaky = 0
    total_hard_fail = 0

    with open(artifact, "w") as fh:
        for case in CASES:
            runs_failed = 0
            first_fail_msgs: List[str] = []
            first_soft: List[str] = []
            fatal_exc: Optional[Exception] = None
            for n in range(repeat):
                try:
                    ctx = build_context(case, rcsb)
                    result = run_assistant(case["query"], ctx)
                    hard = evaluate(result, case.get("expect", {})) + no_tool_leak(result)
                    soft = evaluate(result, case.get("soft", {}))
                    record = {
                        "case": case["id"], "run": n, "query": case["query"],
                        "kind": result.kind, "hard_fail": hard, "soft_warn": soft,
                        "trace": [t.tool for t in result.trace],
                        "result": result.model_dump(),
                    }
                except Exception as e:  # noqa: BLE001
                    hard = [f"exception: {type(e).__name__}: {e}"]
                    soft = []
                    if _is_fatal_provider_error(e):
                        fatal_exc = e
                    record = {
                        "case": case["id"], "run": n, "query": case["query"],
                        "kind": "EXCEPTION", "hard_fail": hard, "soft_warn": soft,
                        "traceback": traceback.format_exc(),
                    }
                fh.write(json.dumps(record, default=str) + "\n")
                if fatal_exc is not None:
                    print("\n!! ABORTED — provider auth/quota error, not an assistant failure:")
                    print(f"   {type(fatal_exc).__name__}: {str(fatal_exc)[:280]}")
                    print(f"   {case_pass} case(s) fully passed before this. artifact {artifact}")
                    return 3
                if hard:
                    runs_failed += 1
                    if not first_fail_msgs:
                        first_fail_msgs = hard
                if n == 0:
                    first_soft = soft

            if runs_failed == 0:
                status = "PASS"
                case_pass += 1
            elif runs_failed == repeat:
                status = "FAIL"
                total_hard_fail += 1
            else:
                status = f"FLAKY {repeat - runs_failed}/{repeat}"
                case_flaky += 1
                total_hard_fail += 1

            print(f"  [{status:>10}] {case['id']}")
            if first_fail_msgs:
                for m in first_fail_msgs:
                    print(f"               - {m}")
            for m in (first_soft or []):
                print(f"               ~ (soft) {m}")

    print(f"\n=== {case_pass}/{len(CASES)} passed"
          + (f", {case_flaky} flaky" if case_flaky else "")
          + f" — artifact {artifact} ===")
    return 1 if total_hard_fail else 0


if __name__ == "__main__":
    sys.exit(main())
