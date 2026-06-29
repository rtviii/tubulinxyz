# MAP partner-protein interface pipeline (Tier 2)

Built 2026-06-18. This documents the MAP-on-tubulin interface pipeline and, for every
layer, how it mirrors or diverges from the existing ligand binding-site pipeline. It
is written for a planned rework: read the "Divergences" and "Things to revisit"
sections first if you are about to change this.

## Mental model

For a ligand, a "binding site" is the protein residues near the ligand in a structure,
aggregated across structures to master-alignment positions and transposed onto any
target. For a MAP (microtubule-associated protein), the "binding site on tubulin" is
the tubulin residues that contact the MAP, aggregated the same way. The MAP plays the
role the ligand plays; the contacted tubulin chain plays the role the protein played.
Everything downstream keys on the master_index of tubulin residues, so the aggregation,
the frequency heatmap, and the frontend transpose are reused unchanged.

The whole point: a MAP interface and a ligand pocket are the same shape of data, so the
read/aggregate/ground/paint half of the stack is shared. Only the data SOURCE differs
(a partner polymer chain instead of a non-polymer ligand), which is the extraction +
storage half.

## Layer-by-layer comparison

| Concern | Ligand pipeline | Partner (MAP) pipeline |
| --- | --- | --- |
| Extraction seed | each non-polymer ligand instance | each protein chain |
| Extraction output | `ligand_neighborhoods` (protein neighbors of a ligand) | `partner_contacts` (protein neighbors of a chain, family-blind, both directions) |
| Classification | n/a (ligand is non-polymer by definition) | collector keeps only MAP-source / tubulin-target pairs |
| Residue payload | `BindingSiteResidue` (tubulin residues) | `BindingSiteResidue` (tubulin residues) -- reused verbatim |
| master_index stamping | `augment_binding_sites` | `augment_partner_contacts` (near-identical) |
| Graph edge | `(NonpolymerInstance)-[:NEAR_POLYMER]->(PolypeptideInstance)` | `(PolypeptideInstance{map})-[:PARTNER_NEAR_POLYMER]->(PolypeptideInstance{tubulin})` |
| Edge payload | `residues_json` (tubulin residues) + `residue_count` | identical |
| Canonical aggregation | `get_canonical_binding_site(chemical_id, family)` | `get_canonical_partner_site(map_family, tubulin_family)` -- keyed by BOTH families |
| Per-structure dedup | per-record (has a frequency>1.0 bug) | per-structure (fixed) |
| Min-count gate | none | `MIN_PARTNER_STRUCTURES = 5` |
| REST | `GET /ligands/canonical-site/{chemical_id}/{family}` | `GET /partners/canonical-site/{map_family}/{tubulin_family}` |
| Per-structure REST | `GET /ligands/neighborhoods/...` (exists) | not built (deferred) |
| Retrieval tool | `get_binding_site` | `get_partner_binding_site` |
| Tool result join key | `family` = the protein family | `family` = the TUBULIN family (map family in separate `map_family` key) |
| Harvest helper | `_binding_chem` | `_interface_partner` |
| Entity category | `binding` (amber) | `interface` (emerald) |
| residue_set key | `chemical_id` | `label` (prettified map family); no chemical_id |
| Expert paint hook | `useCanonicalBindingSite` | not built (deferred) |

## Decisions, by layer

### 1. Extraction (`scripts_and_artifacts/extract_structure_data.tsx`)

`extractPartnerContacts` runs `includeSurroundings({radius: 5, 'as-whole-residues': true})`
around each protein chain (seeds restricted to the observed-sequence chains, so ligand /
water "chains" are never seeded), then buckets the neighbouring protein residues by their
chain. It emits one `PartnerContact{partner_auth_asym_id, contacted_auth_asym_id,
contact_residues}` per (seed chain -> neighbour chain) pair, BOTH directions.

- Decision: simple 5A shell, not typed `InteractionsProvider` (`extract_ixs.tsx`). The
  canonical footprint is a residue set; interaction types (H-bond/hydrophobic/...) would
  be thrown away, and `extract_ixs.tsx` needs the heavy headless gl stack and is single-
  target per run. This is a true 1:1 mirror of `extractLigandNeighborhood`. (Revisit if
  shell noise -- "near but not contacting" -- proves too high; see
  `notes`/`reference_molstar_interchain_bonds.md` for the typed-contacts upgrade path.)
- Decision: extraction is FAMILY-BLIND. The `.tsx` has no HMM/family knowledge (that lives
  only in Python). So it emits all inter-chain protein contacts and the collector decides
  which side is a MAP and which is tubulin. Both directions are emitted; the collector
  keeps only MAP->tubulin.
- Radius 5 + as-whole-residues matches the ligand pass exactly.

### 2. Types (`lib/types.py`)

- `BindingSiteResidue` is REUSED for partner residues (it already carries `master_index`
  and `to_dict()`, and `populate_by_name=True` lets us construct by field name).
- Two new types, because classification happens between them:
  - `RawPartnerContact` = one directed family-blind pair from extraction
    (`partner_auth_asym_id`, `contacted_auth_asym_id`, `residues` on the contacted chain).
    Held on `MolstarExtractionResult.partner_contacts`.
  - `PartnerContact` = the classified, per-MAP-chain result (`partner_auth_asym_id` = the
    MAP chain, `partner_family` = the MAP family enum string, `residues` = tubulin-side
    residues with master_index). Held on `TubulinStructure.partner_contacts`. This is the
    analogue of `LigandBindingSite`.
- `PartnerContactsFile` mirrors `LigandBindingSitesFile`.

### 3. Classification + augmentation (`collector.py`, `augmentation.py`)

- `TubulinETLCollector._build_partner_contacts` filters the raw contacts to the direction
  where the seed is a MAP (`is_map_family`) and the contacted chain is a tubulin chain
  present in `chain_mappings` (`is_tubulin_family` + master-mapped). It unions the tubulin
  residues per MAP chain into one `PartnerContact`. So a MAP that touches alpha and beta
  produces ONE `PartnerContact` whose residues span both tubulin chains.
- Decision: tubulin-side only. MAP-side residues are never collected (no master alignment
  exists for MAP chains; `chain_mappings` holds tubulin chains only). This is what we paint
  on a tubulin target.
- `augment_partner_contacts` is a near-verbatim copy of `augment_binding_sites` -- it
  stamps `master_index` on the tubulin residues via
  `chain_mappings[auth_asym_id].auth_seq_id_to_master`. Same maps the ligand path uses.

### 4. Storage (`neo4j_tubxz/node_partner_site.py`, `db_lib_builder.py`)

- New edge `(:PolypeptideInstance{map})-[:PARTNER_NEAR_POLYMER {residues_json,
  residue_count}]->(:PolypeptideInstance{tubulin})`. One edge per (MAP chain, tubulin
  chain); `residues_json` is the tubulin residues for that chain via `to_dict()`.
- Decision: a DISTINCT edge label, never `NEAR_POLYMER`. The ligand queries
  (`get_canonical_binding_site`, `get_ligand_neighborhoods_for_polymer`,
  `PolypeptideEntityQueryBuilder.has_ligand_ids`) assume the left endpoint is a
  `NonpolymerInstance` backed by a `Chemical`; a MAP edge on that label would corrupt them.
- Decision: directed MAP -> tubulin, tubulin-side payload only.
- `MERGE` -> idempotent; re-running overwrites `residues_json`/`residue_count`.
- Ingested in `add_total_structure` Phase 5b, mirroring the ligand Phase 5.

### 5. Canonical aggregation (`db_lib_reader.py`)

`get_canonical_partner_site(map_family, tubulin_family)` copies
`get_canonical_binding_site` with two deliberate changes:

- Keyed by BOTH families: `$map_family` on the source entity, `$tubulin_family` on the
  contacted entity. Necessary because a MAP contacts alpha and beta, and alpha/beta have
  SEPARATE master alignments -- a result must be scoped to the tubulin numbering being
  painted. (The ligand reader keys by `chemical_id` + one family.)
- Per-structure dedup FIX. The ligand reader counts per-record (per edge), so when one
  structure contacts the same master_index through multiple instances/chains, the count
  exceeds the structure count and frequency can exceed 1.0. MAPs hit that case far more
  often (oligomeric, multiple copies, contacting two tubulin chains). This copy groups
  master indices by `rcsb_id` first, so each master_index is counted once per structure
  and frequency stays in [0, 1]. Verified on the stathmin smoke set: max frequency = 1.0.
- Gate: returns `None` below `MIN_PARTNER_STRUCTURES = 5`. A single observation is not a
  "canonical" site. Single source of truth -- both the REST 404 and the tool's
  `found:false` inherit it. The threshold is a guess; tune it.
- Reuses `CanonicalBindingSiteResidue` (identical `{master_index, count, frequency}`), so
  the frontend heatmap path is untouched. New wrapper model `CanonicalPartnerSite`.

### 6. REST (`api/routers/router_partners.py`, `main.py`)

`GET /partners/canonical-site/{map_family}/{tubulin_family}`, `operation_id =
get_partner_canonical_site`, 404 below the gate / when no data. Registered with the
`/partners` prefix via the direct-import convention (like `router_annotations`). Path
params pass straight through; families are NOT uppercased (map_*/tubulin_* are lowercase
enums; only chemical_id/rcsb_id get uppercased in the data layer).

Decision: no per-structure neighborhoods endpoint (the ligand path has one). The canonical
endpoint serves the landing-assistant footprint; per-structure painting is an expert-mode
feature, deferred.

### 7. Retrieval tool + harvest (`retrieval.py`, `orchestrator.py`)

- `get_partner_binding_site(map_family, tubulin_family)` mirrors `get_binding_site`.
  CRITICAL: the result's `family` key is the TUBULIN family, because the harvest joins
  `result['family']` against the demo's chain families (`fam_to_chain`). The MAP family
  lives in a separate `map_family` key. (Putting `map_*` in `family` would silently emit
  nothing -- no demo chain has a map family.)
- Harvest: `get_partner_binding_site` added to `_HARVEST_TOOLS`, `_TOOL_CATEGORY` maps it
  to `interface`. A new `_interface_partner` helper (parallel to `_binding_chem`) drives a
  `residue_set` emission branch keyed by `map_family` with `category="interface"` and
  `label` = prettified map family (no `chemical_id`). The `residue_range`/`region`
  harvest blocks already tag with `category` from `_TOOL_CATEGORY`, so they come out
  `interface` for free.
- Harvest runs ONLY on the landing page (`ctx.page == 'landing'`). On other pages the
  model emits actions itself; there is no `interface_contacts` AddAnnotationTrack kind yet
  (would be needed for structure/catalogue-page auto-painting).
- `EntityRef.category` is free-text, so `interface` needed no schema change.

### 8. Frontend (`fend_tubulinxyz`)

- `interface` added to TWO color maps: `entityHighlight.ts` `CATEGORY_TONE` (emerald
  pill tint) and `src/lib/colors/annotationPalette.ts` `CATEGORY_PAINT` (`#10b981` 3D
  label). Both are needed; adding only one fails silently. Emerald chosen to be distinct
  from binding/amber, modification/indigo, variant/orange, and to avoid retired teal.
- `applyHighlight` and the `EntityRef` type needed NO change -- a partner interface is just
  `EntityRef{kind:'residue_set', auth_asym_id, positions, category:'interface'}`.
- No `yarn codegen` needed for this scope: the landing assistant paints via harvested
  entities the frontend already renders; nothing on the client calls `/partners` directly.
  (A `yarn codegen` would add a `getPartnerCanonicalSite` hook for a future expert hook.)

### 8b. Prompt + eval

- The "no MAP binding-site tool yet" prelude line was replaced with guidance to call
  `get_partner_binding_site`, and a READ TOOLS prose line was added.
- `scripts/eval_assistant.py`: `landing_demo_honesty` category whitelist gained
  `interface`; `eb1_binds_where_honest` was upgraded to `eb1_binds_where_interface`
  (now expects the tool call); a landing positive case `stathmin_interface_landing` was
  added. The eval makes ~18 LLM calls and was not run in this session; run it post-backfill.

## Divergences from the ligand pipeline (the non-mirror parts)

1. Net-new extraction. The ligand path already extracted neighborhoods; the partner path
   added a whole new family-blind protein-protein pass.
2. A classification step (MAP vs tubulin) sits between extraction and storage; ligands have
   none.
3. Distinct edge label (`PARTNER_NEAR_POLYMER`) to protect the ligand queries.
4. Canonical query keyed by TWO families, not chemical_id + one.
5. Per-structure dedup is FIXED here but the ligand reader still has the per-record bug.
6. A min-structure gate exists here; ligands have none.
7. Tool result `family` is the tubulin (contacted) family, an alias chosen for the harvest
   join; the ligand tool's `family` is the contacted family too, but there it is also the
   semantic family -- here the semantic key (the MAP) is separate.

## Things to revisit in the rework

- Simple 5A shell vs typed interactions. If the shell admits too many "near but not
  contacting" residues, switch the extraction to `InteractionsProvider` (see
  `extract_ixs.tsx`, currently experimental/unused) and/or filter by contact count.
- Fix the per-record dedup bug in the LIGAND `get_canonical_binding_site` too (we only
  fixed the partner copy). Same per-rcsb_id grouping applies.
- `MIN_PARTNER_STRUCTURES = 5` is a guess. Tune per family; some well-populated families
  could justify a lower bar, sparse ones a higher one.
- Per-structure partner neighborhoods endpoint + a `usePartnerBindingSite` expert hook are
  deferred -- add them when wiring expert-mode painting. The expert hook can mirror
  `useCanonicalBindingSite` almost verbatim (only the endpoint URL + base color differ);
  note `applyColorscheme` calls `restoreDefaultColors()` internally, so co-displaying a
  partner footprint and a ligand pocket needs a merged colorings array, not two calls.
- Symmetry: MAP-side residues are not stored. If a MAP-centric view is wanted later, add a
  second edge/property; do not overload the tubulin-side edge.
- Structure-page auto-painting: there is no `interface_contacts` AddAnnotationTrack kind,
  so the interface only auto-appears on the landing demo. Add a track kind for the
  structure/catalogue pages if wanted.
- Disk profiles: the drift-free backfill (below) writes edges to the graph but does NOT
  write `partner_contacts` to the profile JSON. The canonical endpoint reads the graph, so
  the feature works, but the disk profiles are not self-describing. A real re-collect
  (`collect-one -o`) does write them.

## Backfill (how the edges get populated)

The feature is inert until `PARTNER_NEAR_POLYMER` edges exist. Two paths:

- Full re-collect: `cli.py collect-all -o` + `upload-all -o`. Simplest (no new code) and
  also refreshes profiles to current backend logic, BUT re-fetches RCSB, so it incidentally
  changes non-Tier-2 fields. Observed on 9MLF: `entities[*].uniprot_accessions` emptied
  (`["A0A287AGU7"]`->`[]`), and `nonpolymers`/`assembly_map` element ordering churned. Not
  purely additive.
- Drift-free additive (`scripts/backfill_partner_contacts.py`, used here): extracts to a
  temp file, rebuilds the MAP/tubulin classification + augmentation inputs from the EXISTING
  profile (its families + `chain_index_mappings`), and adds ONLY the edges via
  `process_all_partner_sites`. No RCSB refetch, no profile rewrite, additive + reversible
  (`--delete-all`). Reuses the production `_build_partner_contacts` / `augment_partner_contacts`
  so the edges match a real re-collect. Usage:

  ```
  set -a; . ./.env; set +a
  ./venv/bin/python scripts/backfill_partner_contacts.py --family map_stathmin --limit 8
  ./venv/bin/python scripts/backfill_partner_contacts.py --maps-only
  ./venv/bin/python scripts/backfill_partner_contacts.py --delete-all
  ```

## Backfill results (2026-06-18, drift-free additive path, 7 workers)

Populated 11,306 `PARTNER_NEAR_POLYMER` edges across 594 structures (of the 704
has_any_map set; the other ~110 yield nothing -- see EB note). Parallel run over the
578 not-already-done took ~7 minutes; 0 errors.

Canonical sites now available (structure_count, by tubulin family):
- map_stathmin: beta 332 / alpha 332 (full coverage)
- map_ttll_glutamylase_short: beta 292 / alpha 292 (full)
- map_kinesin13: beta 120 / alpha 132
- map_gcp2_3 / map_gcp4 / map_gcp5_6: against tubulin_GAMMA (40 / 28 / 29) -- NOT
  alpha/beta. GCPs are gamma-TuRC; always query them with tubulin_gamma.
- map_doublecortin: ~12; map_numa: 12; map_rib72_efhc / map_fap20_cfap20 / map_pacrg /
  map_cfap53 / map_nme7: 17-25; map_tau: 6; map_clip170: ~6.
- Max frequency stays <= 1.0 everywhere (per-structure dedup verified). Enzymes /
  transient binders (TTLL, clip170, spag6) show max_freq < 1.0, as expected.

Two findings for the rework:
- GCP -> tubulin_gamma: the assistant prompt now hints which tubulin family each MAP
  contacts (GCPs -> gamma; EB/stathmin/kinesin/etc -> beta/alpha). Without the hint the
  model would default to beta and get an honest "below gate".
- EB coverage gap: only 8 of 103 map_eb_family structures have a master-mapped tubulin
  entity in their EXISTING profile; the other 95 have unclassified / unmapped tubulin
  (family=None) or are EB-only. So EB's canonical site rests on 8 structures. This is an
  upstream profile classification/alignment coverage issue, NOT Tier 2 -- a real
  re-collect with the current HMM/aligner may recover many. Worth checking during the
  rework whether other big families (some kinesin13/eb microtubule sets) are similarly
  under-counted because their lattice tubulin didn't classify.

## Verify recipes

- Canonical query (no LLM): `db_reader.get_canonical_partner_site('map_stathmin', 'tubulin_beta')`
  -> structure_count + residues with frequency in [0,1].
- REST: `curl /partners/canonical-site/map_stathmin/tubulin_beta`.
- Tool: `run_retrieval_tool('get_partner_binding_site', {'map_family':..., 'tubulin_family':...})`
  -> check `result['family']` is the tubulin family.
- End-to-end (LLM): `run_assistant('Where does stathmin contact tubulin?', LANDING_DEMO)`
  -> expect `get_partner_binding_site` in the trace and `residue_set`/`residue_range`
  entities with `category='interface'` grounded on demo chains A/B.
- Eval: `python scripts/eval_assistant.py` (~18 LLM calls; needs Neo4j + the edges).
