#!/usr/bin/env python
"""
scripts/backfill_partner_contacts.py

DRIFT-FREE, additive backfill of MAP->tubulin PARTNER_NEAR_POLYMER edges (Tier 2).

WHY THIS EXISTS / HOW IT DIFFERS FROM A FULL RE-COLLECT
-------------------------------------------------------
A full `cli.py collect-all -o` + `upload-all -o` would populate the same edges,
but it re-fetches every structure from the RCSB GraphQL API and rewrites every
profile JSON. That re-fetch incidentally changes non-Tier-2 content (observed on
9MLF: `entities[*].uniprot_accessions` emptied `["A0A287AGU7"]`->`[]`, and
`nonpolymers` / `assembly_map` element ordering churned). None of that is related
to the partner-interface feature.

This script avoids the drift entirely:
  1. Runs the Molstar extraction to a TEMP file (the CIF is already on disk; no
     network). This yields the raw, family-blind inter-chain `partner_contacts`.
  2. Rebuilds the MAP/tubulin classification + master-index augmentation inputs
     from the EXISTING on-disk profile (its `entities[*].family` and
     `entities[*].chain_index_mappings`) -- it does NOT re-classify or re-align.
  3. Adds ONLY the PARTNER_NEAR_POLYMER edges to the graph via
     `process_all_partner_sites` (MERGE -> idempotent). No profile is rewritten;
     RCSB is never queried. The single graph mutation is additive + reversible.

It reuses the SAME production functions the real collector uses
(`TubulinETLCollector._build_partner_contacts`, `augment_partner_contacts`,
`process_all_partner_sites`), so the edges are identical to what a real
re-collect+upload would produce -- minus the drift.

CAVEATS (documented for the planned rework):
  - The disk profiles do NOT get `partner_contacts` written here. The canonical
    endpoint + assistant tool read edges from the GRAPH (db_reader), not the
    profile JSON, so the feature works without it. A future "real" backfill /
    rework can persist `partner_contacts` to disk too (run collect-one -o, which
    already does, once the RCSB drift is acceptable or pinned).
  - This uses the existing profile's classification + alignment as-is. If a
    profile predates a classification/alignment improvement, the edges reflect
    the old mapping -- but that is exactly consistent with the served data.
  - Idempotent: MERGE overwrites residues_json/residue_count on re-run.

USAGE
-----
  set -a; . ./.env; set +a
  ./venv/bin/python scripts/backfill_partner_contacts.py --family map_stathmin --limit 8
  ./venv/bin/python scripts/backfill_partner_contacts.py --maps-only
  ./venv/bin/python scripts/backfill_partner_contacts.py --ids 9MLF 8IBN --dry-run
  ./venv/bin/python scripts/backfill_partner_contacts.py --delete-all   # remove all partner edges
"""
import argparse
import os
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# `python scripts/backfill_partner_contacts.py` sets sys.path[0] to scripts/, not
# the repo root -- add the root so `api` / `lib` / `neo4j_tubxz` import.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from loguru import logger

from api.config import PROJECT_ROOT
from lib.etl.assets import TubulinStructureAssets
from lib.etl.collector import TubulinETLCollector
from lib.etl.molstar_bridge import run_molstar_extraction
from lib.etl.augmentation import augment_partner_contacts
from lib.etl.classification import is_tubulin_family
from lib.types import PolypeptideEntity
from neo4j_tubxz.db_lib_reader import db_reader
from neo4j_tubxz.models import StructureFilters
from neo4j_tubxz.node_partner_site import process_all_partner_sites

EXTRACT_SCRIPT = PROJECT_ROOT / "scripts_and_artifacts" / "extract_structure_data.tsx"


def _all_ids(filters: StructureFilters) -> list[str]:
    """Paginate list_structures to collect every matching rcsb_id."""
    ids: list[str] = []
    cursor = None
    while True:
        f = filters.model_copy(update={"cursor": cursor, "limit": 200})
        resp = db_reader.list_structures(f)
        ids.extend(s.rcsb_id for s in resp.data)
        if not resp.has_more or not resp.next_cursor:
            break
        cursor = resp.next_cursor
    return ids


def resolve_targets(args) -> list[str]:
    if args.ids:
        return [i.upper() for i in args.ids]
    if args.family:
        return _all_ids(StructureFilters(has_polymer_family=[args.family]))
    if args.maps_only:
        return _all_ids(StructureFilters(has_any_map=True))
    raise SystemExit("Specify one of --ids, --family, or --maps-only (or --delete-all).")


def backfill_one(rcsb_id: str, tmp_dir: Path, dry_run: bool) -> dict:
    """Extract -> classify+augment from existing profile -> add edges. Returns a
    small result dict; never raises (errors are captured per-structure)."""
    rcsb_id = rcsb_id.upper()
    out = {"rcsb_id": rcsb_id, "contacts": 0, "edges": 0, "status": "ok"}
    try:
        collector = TubulinETLCollector(rcsb_id)
        cif = Path(collector.assets.paths.cif)
        if not cif.exists():
            out["status"] = "skip:no_cif"
            return out
        if not os.path.exists(collector.assets.paths.profile):
            out["status"] = "skip:no_profile"
            return out

        profile = collector.assets.profile()

        # Extraction to a temp path (CIF is local; no network, no profile rewrite).
        temp_out = tmp_dir / f"{rcsb_id}_partner_raw.json"
        mres = run_molstar_extraction(
            cif_path=cif, rcsb_id=rcsb_id, output_path=temp_out,
            script_path=EXTRACT_SCRIPT, project_root=PROJECT_ROOT,
        )
        if not mres:
            out["status"] = "skip:extraction_failed"
            return out

        # Rebuild classification + augmentation inputs from the EXISTING profile.
        chain_to_entity = {p.auth_asym_id: p.entity_id for p in profile.polypeptides}
        entity_families = {
            eid: ent.family
            for eid, ent in profile.entities.items()
            if isinstance(ent, PolypeptideEntity)
        }
        # chain_mappings: tubulin chains only (mirrors collector Phase 3). The
        # stored ChainIndexMappingData exposes .auth_seq_id_to_master, which is all
        # _build_partner_contacts (membership) and augment_partner_contacts (lookup)
        # need -- duck-typed identically to the live ChainIndexMapping.
        chain_mappings = {}
        for ent in profile.entities.values():
            if isinstance(ent, PolypeptideEntity) and is_tubulin_family(ent.family):
                for aid, cmd in (ent.chain_index_mappings or {}).items():
                    chain_mappings[aid] = cmd

        # Reuse the production classification + augmentation, verbatim.
        contacts = collector._build_partner_contacts(
            raw_contacts=mres.partner_contacts,
            chain_to_entity=chain_to_entity,
            entity_families=entity_families,
            chain_mappings=chain_mappings,
        )
        contacts = augment_partner_contacts(contacts, chain_mappings)
        out["contacts"] = len(contacts)

        if dry_run:
            out["status"] = "dry-run"
            return out

        with db_reader.adapter.driver.session() as s:
            out["edges"] = s.execute_write(process_all_partner_sites(rcsb_id, contacts))
    except Exception as e:  # per-structure resilience
        out["status"] = f"error:{type(e).__name__}:{e}"
    return out


def delete_all_partner_edges() -> int:
    with db_reader.adapter.driver.session() as s:
        rec = s.run(
            "MATCH ()-[r:PARTNER_NEAR_POLYMER]->() WITH r LIMIT 1000000 "
            "DELETE r RETURN count(r) AS c"
        ).single()
        return rec["c"] if rec else 0


def existing_partner_structures() -> set:
    """rcsb_ids that already have at least one PARTNER_NEAR_POLYMER edge."""
    with db_reader.adapter.driver.session() as s:
        return {
            r["rid"]
            for r in s.run(
                "MATCH (mi:PolypeptideInstance)-[:PARTNER_NEAR_POLYMER]->() "
                "RETURN DISTINCT mi.parent_rcsb_id AS rid"
            )
        }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--ids", nargs="+", help="Explicit rcsb_ids")
    g.add_argument("--family", help="A MAP family enum, e.g. map_stathmin")
    g.add_argument("--maps-only", action="store_true", help="All structures with any MAP")
    ap.add_argument("--limit", type=int, default=None, help="Cap the number of structures")
    ap.add_argument("--workers", type=int, default=1, help="Concurrent extraction workers (threads; each shells one Molstar subprocess)")
    ap.add_argument("--skip-existing", action="store_true", help="Skip structures that already have PARTNER_NEAR_POLYMER edges")
    ap.add_argument("--dry-run", action="store_true", help="Extract + classify but write no edges")
    ap.add_argument("--delete-all", action="store_true", help="Delete ALL PARTNER_NEAR_POLYMER edges and exit")
    args = ap.parse_args()

    if args.delete_all:
        n = delete_all_partner_edges()
        logger.info(f"Deleted {n} PARTNER_NEAR_POLYMER edges.")
        return

    targets = resolve_targets(args)
    if args.skip_existing:
        done = existing_partner_structures()
        before = len(targets)
        targets = [t for t in targets if t not in done]
        logger.info(f"--skip-existing: {before - len(targets)} of {before} already have edges, skipping.")
    if args.limit:
        targets = targets[: args.limit]
    n = len(targets)
    logger.info(f"Backfilling partner contacts for {n} structures (workers={args.workers}, dry_run={args.dry_run}).")

    totals = {"ok": 0, "edges": 0, "contacts": 0, "skipped": 0, "errors": 0}

    def record(i: int, r: dict) -> None:
        st = r["status"]
        if st in ("ok", "dry-run"):
            totals["ok"] += 1
            totals["edges"] += r["edges"]
            totals["contacts"] += r["contacts"]
        elif st.startswith("skip"):
            totals["skipped"] += 1
        else:
            totals["errors"] += 1
        logger.info(f"[{i}/{n}] {r['rcsb_id']}: {st} (contacts={r['contacts']}, edges={r['edges']})")

    with tempfile.TemporaryDirectory(prefix="partner_backfill_") as td:
        tmp_dir = Path(td)
        if args.workers > 1:
            # Each thread shells its own Molstar extraction subprocess and opens its
            # own Neo4j session; distinct structures touch distinct graph nodes, so
            # there is no cross-structure write contention. Mirrors render_thumbnails.
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                futs = [ex.submit(backfill_one, rid, tmp_dir, args.dry_run) for rid in targets]
                for i, fut in enumerate(as_completed(futs), 1):
                    record(i, fut.result())
        else:
            for i, rid in enumerate(targets, 1):
                record(i, backfill_one(rid, tmp_dir, args.dry_run))

    logger.info(
        f"DONE. ok={totals['ok']} skipped={totals['skipped']} errors={totals['errors']} "
        f"| total partner contacts={totals['contacts']} edges={totals['edges']}"
    )


if __name__ == "__main__":
    main()
