Summary for your meeting
Honest tier-based assessment:

Tier 1 - production quality after the terminal deletion fix: alpha, beta. The core logic is sound, the MSAs are good, the remaining artifacts are localized to the tails and are fixable with the reclassification approach we discussed.

Tier 2 - gamma is partially working. The genuine human gamma-1 and gamma-2 from gamma-TuRC structures are being found and aligned, but the annotations are polluted by (a) 19 systematic internal spurious deletions from MSA column bias and (b) misclassified co-complex proteins. Needs MSA column filtering and HMM score thresholding to fix.

Tier 3 - delta and epsilon are not usable yet. The MSAs need to be rebuilt with broader phylogenetic coverage before these families can be annotated reliably. Flag this as known future work.


# Isotype-aware substitution calling via UniProt references

# Problem being solved

The current substitution calling compares every aligned residue against the MSA consensus - a single character per position computed by plurality vote across all members of the family alignment. This means positions where isotypes genuinely differ from each other get flagged as substitutions. For example, if 4 out of 6 beta isotypes have S at position 48 but TUBB1 consistently has N there, every TUBB1 structure gets a spurious S->N substitution call. This is noise that obscures real mutations.

# Core idea

Keep the MSA alignment and master_index mapping exactly as-is - that's the cross-structure coordinate system and it works. Change only what wild_type means when emitting a substitution: instead of the consensus character, use the residue at the corresponding position in the UniProt canonical sequence for that specific isotype.
Implementation sketch
A UniProtReferenceLibrary class initialized from the audit JSON, holding reviewed sequences keyed by accession and family. It exposes two methods: get_reference(accession) for direct lookup, and find_closest(canonical_sequence, family) for inference via pairwise identity. The library is built once at startup. In sequence_alignment.py, before calling substitutions, the aligner looks up the appropriate reference sequence and uses it instead of self.consensus.get_residue() at each aligned position.

# Edge case cascade

- Reviewed UniProt accession present on the entity: use it directly. Full confidence. This covers 99% of your current dataset.
- Unreviewed accession present: use its sequence as reference but flag the substitution calls as lower-confidence. Optionally verify the isotype assignment by pairwise identity against the reviewed references.
- No accession (6 entities currently): run pairwise identity against all reviewed references for that family, assign to the closest, flag as inferred. Since most of these are from organisms without many paralogous isotypes, ambiguity is unlikely in practice.
- Accession maps to a non-tubulin (EFHC2, TUBGCP2 etc.): these indicate misclassification upstream - the entity shouldn't have been processed through alignment at all. These should be caught by a post-classification filter before alignment runs.


The confidence tier of the isotype assignment propagates to the source field on each SequenceVariant, so downstream consumers know which substitution calls to trust.

What this buys you for the presentation

The substitution counts drop dramatically for well-studied isotypes - the 20-30 substitutions you currently see on many beta structures will reduce to only genuine within-isotype deviations. The calls that remain will be directly comparable to literature mutation databases including Morisette, because both are now anchored to the same per-isotype UniProt sequence.
-------

## Pipeline audit summary - tubulin ETL system

### Architecture overview

The pipeline processes PDB structures of tubulin through five phases: raw data acquisition (RCSB API + Molstar extraction), HMM-based family classification, sequence alignment against family MSAs, binding site augmentation, and profile assembly.

The central design choice is a two-sequence architecture:
- **canonical sequence**: `one_letter_code_can` from the RCSB entity, sourced from SEQRES records. Represents what the depositor claims was in the construct, including residues with no electron density. Used for alignment and variant calling.
- **observed sequence**: extracted by Molstar from ATOM records only. Carries both `auth_seq_id` and `label_seq_id` per residue. Used for classification and for building per-chain coordinate mappings.

The bridge between the two is `label_seq_id`, which by mmCIF convention is the 1-based index into the SEQRES sequence. Molstar reads this directly from `_atom_site.label_seq_id`. This has been verified correct from the extraction script and confirmed by the 2WBE sample where the gap at positions 35-60 appears in the observed sequence but not in the canonical, and correctly produces `unresolved_positions` rather than deletion calls.

The master alignment (MSA) per family serves as the coordinate lingua franca. Every canonical sequence is profile-aligned against the family MSA using MUSCLE, producing a `label_seq_id -> master_index` mapping. Per-chain `auth_seq_id -> master_index` mappings are then derived cheaply from that entity-level result plus the observed residue data. All downstream annotations (variants, binding sites, PTMs, interfaces) are stored with `master_index` so any two structures in the same family can be compared positionally.

---

### What is sound

- The canonical/observed split correctly separates genetic variants from unresolved density. Unresolved loops and disordered tails that are present in the canonical sequence but absent from the observed sequence land in `unresolved_positions`, not in the deletion list. This is the core architectural claim and it holds.
- The `label_seq_id` bridge is reliable for standard structures. Molstar reads the mmCIF field directly and the extraction script has been verified.
- HMM classification is working well for alpha and beta. Scores are strong, the family assignments are reliable.
- UniProt accession coverage is 99% across all tubulin entities (1529/1535), with all 102 unique accessions having sequences successfully retrieved. Alpha and beta are dominated by reviewed Swiss-Prot entries.
- Per-entity invariance of the canonical sequence is guaranteed by PDB: if any chain instance of an entity differs in sequence, it gets its own entity. So using one canonical sequence per entity is correct.

---

### Confirmed bugs

**X residue substitutions.** `one_letter_code_can` can contain `X` for ambiguous residues. The current code emits a substitution `wildtype->X` at those positions. These are meaningless and should be filtered out before variant calling. Fix: skip any column where `t_char == 'X'` in `_build_entity_result`.

**Insertion code collision in Molstar extraction.** The deduplication key in `extract_structure_data.tsx` is `${authAsymId}:${authSeqId}`. For residues with insertion codes (100A, 100B), Molstar returns only the integer part of `auth_seq_id`, so two distinct residues can silently collapse to one. Rare in tubulin but worth knowing about.

---

### Deletion noise - diagnosis and fix

**Mechanism.** Terminal truncations in PDB constructs - commonly the disordered C-terminal tails of tubulin - produce alignment gaps at the ends of the canonical sequence. Because the MSA consensus has non-gap characters at those positions (contributed by species with longer tails), the code emits deletion variants there. These are technically correct descriptions of the construct but are biologically misleading and constitute most of the deletion noise.

**Measured impact from diagnostic script across full dataset:**
- Beta: 100% of all deletions are terminal, avg 9.2 per entity. The N-cluster consistently at master positions 431-432, C-cluster varying between 434-452 depending on the isotype tail length. Zero internal deletions worth worrying about.
- Alpha: 89% terminal, avg 2.5. 205 internal deletions across 736 entities - acceptable.
- Gamma: 43% terminal, but 1054 internal deletions across 55 entities. Two distinct problems (see below).
- Delta: 85% internal, avg 73 total. MSA failure, not structure problems.
- Epsilon: 84% internal, avg 128 total. MSA failure.

**Fix for alpha/beta.** After building variants in `_build_entity_result`, detect contiguous deletion runs at the N and C ends of the canonical coverage range and reclassify them as `VariantType.TERMINAL_TRUNCATION`. Logic: find the lowest and highest canonical positions with a non-None master_index mapping, find their corresponding master positions, then any deletion whose master_index falls outside that range is terminal. Requires adding `TERMINAL_TRUNCATION` to `VariantType` enum and optionally an `is_terminal` flag on `SequenceVariant`.

**Gamma-specific problem 1: systematic 19 internal spurious deletions.** All vertebrate gamma structures show exactly 19 internal deletions, consistently at the same master positions. These correspond to MSA columns where non-vertebrate sequences (yeast, worm, plant) contribute residues but all vertebrate sequences have gaps. The consensus calculator puts a non-gap character at those columns, and every vertebrate gamma structure gets 19 deletion calls that are artifacts of the phylogenetically broad MSA. Fix: identify and mask columns where gap frequency among vertebrate sequences exceeds a threshold (e.g. 80%) before computing the consensus, or rebuild the gamma MSA from vertebrate-only sequences.

**Gamma-specific problem 2: misclassified co-complex proteins.** Several structures (9G3X, 9G3Y, 9H9P, 9I8G, 9I8N etc.) have a second entity with 400+ substitutions and 400+ insertions - these are GCP2, GCP3, GCP4, GCP6 subunits of the gamma-TuRC complex being misassigned to `tubulin_gamma` by the HMM. The genuine gamma chain in the same structure is found correctly. Fix: HMM score threshold filtering - the misassigned entities almost certainly have much lower scores than genuine gamma chains.

**Delta and epsilon MSAs are not usable.** Delta MSA has 5 sequences, all vertebrates. Epsilon has 2 sequences, both mammalian. The average substitution counts (309 for delta, 275 for epsilon) compared to alpha (14) confirm that the alignments are producing near-random noise for anything other than the specific organisms in the MSA. Additionally, the delta and epsilon UniProt accessions in the audit revealed a separate misclassification problem: the accessions returned are prokaryotic tubulin homologs (tubZ, cetZ, ftsZ from Bacillus, Methanothrix, Haloferax, Methanopyrus) and phage-encoded proteins, not eukaryotic delta/epsilon tubulin. The HMMs need to be rebuilt with proper eukaryotic sequences.

---

### Substitution calling - the isotype problem and fix

**Current problem.** The MSA consensus is a single character per position across all members of the family. Positions where isotypes genuinely differ from each other are flagged as substitutions relative to consensus. For example TUBB1 structures consistently show `S->N` at position 48, `S->A` at 55, `G->N` at 57 etc. - these are not mutations, they are the expected TUBB1 sequence. This inflates substitution counts and obscures real mutations.

**Proposed fix: isotype-aware substitution calling via UniProt references.** The master_index mapping stays unchanged - it's the coordinate system and works correctly. Only the source of `wild_type` changes. Instead of `consensus.get_residue(master_idx)`, look up the residue at that position in the UniProt canonical sequence for the entity's specific isotype.

Implementation requires a `UniProtReferenceLibrary` class initialized from the audit JSON, holding sequences keyed by accession. Before alignment, build a `uniprot_pos_to_master` mapping by aligning the UniProt sequence against the MSA profile (this is computed once per accession and cached). During variant calling, use `uniprot_sequence[uniprot_pos]` as `wild_type` instead of the consensus character.

Edge case cascade:
1. Reviewed UniProt accession on entity: use directly. Full confidence. Covers 99% of current dataset.
2. Unreviewed accession: use its sequence, flag variants as lower-confidence. Optionally verify by pairwise identity against reviewed references.
3. No accession (6 entities): pairwise identity against all reviewed references for that family, assign to closest, flag as inferred.
4. Accession maps to non-tubulin (EFHC2, TUBGCP2, tubZ etc.): entity should be filtered out at post-classification before reaching alignment. These are upstream misclassifications.

The confidence tier propagates to the `source` field on `SequenceVariant` so downstream consumers know which calls to trust. Substitutions are still reported as mutations relative to the isotype reference - including differences between the PDB canonical and UniProt sequence that may represent genuine engineered mutations in the construct.

---

### Insertions - not yet analyzed

Insertions have not been examined in detail. From the aggregate stats: alpha has avg 3.4 insertions per entity (2509 total), beta has avg 0.3 (196 total), gamma has avg 80.8 (4445 total) which is clearly pathological and tied to the MSA quality problems. For alpha/beta the insertion calls come from MUSCLE adding new columns to accommodate the target sequence - these are positions in the canonical with no corresponding master position, recorded with `canonical_to_master[pos] = None`. The main open question is whether any of these represent real biological insertions versus alignment artifacts at the ends of the sequence or around the internal gap columns. This needs a similar diagnostic pass to what was done for deletions.

---

### MSA quality tier summary

- Alpha, beta: production quality. Good phylogenetic spread, clean alignments, tails are the only systematic problem.
- Gamma: partially working. Core body alignment is reliable for vertebrate sequences but the tail and the 19-column internal gap problem need fixing.
- Delta: needs MSA rebuild with eukaryotic sequences. Current MSA is vertebrate-only and HMM is misclassifying prokaryotic homologs.
- Epsilon: same as delta, worse - only 2 sequences.

---

### Roadmap items in priority order

1. `X` residue substitution filter - confirmed bug, small fix.
2. Terminal truncation reclassification - immediate impact on alpha/beta output quality.
3. Isotype-aware substitution calling via UniProt references - most impactful correctness improvement for alpha/beta.
4. Post-classification filter to reject non-tubulin entities (TUBGCP2, EFHC2, tubZ etc.) before alignment runs.
5. Gamma MSA column masking for the 19-column vertebrate gap problem.
6. Gamma HMM score threshold to reject GCP misassignments.
7. Insertion diagnostic pass analogous to the deletion analysis done above.
8. Delta and epsilon MSA rebuild with eukaryotic sequences.
9. Connect to Morisette et al. mutation database - translation shim from their coordinate system to master_index. (Next topic.)