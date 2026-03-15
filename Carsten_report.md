
# Carsten meet/Thu the 12th


# Feedback

- Monomer view -- make the species/metadata more prominent (probably in the monome row)

# Whether calling substitutions against the uniprot reference is a good idea.

- Highlight isotype differences (including the reference ) -- work out the tiered colorscheme for dis/similar amino acids (Carsten's black/white/red TUB8 example)
- make it easy to add free text/fasta sequeences to alignw
- seqeuence viewer interaction: display 3d b&s amino acid and a metadata Modal for the given residue (ephemeral)


-----

<!-- - spurious termini "deletions"/whether just dropping them as "truncations" is a good compromise -->
<!-- - Annotations confidence tiers -->
<!-- - Whether we want to keep our MSAs focused on core vertebrate species or cover the whole phylogentic tree like Morisette does -->

- Quote for Pascale:

  send the invoice as a proforma/quote
  order number on the invoice --> put on the invoice after agreement



## Context

We want to be able to say: the lysine at this functionally important position in structure A is the same position as the arginine in structure B, and here is the full history of what's been seen at that position across all structures in the database. That requires a shared coordinate system. 


## What we get from the PDB, and why it's not enough

## Family classification

- alpha/beta works well
- gamma works somewhat but produces some false positives (should be fixed by making the MSA more representative)

##  The master alignment as lingua franca

- works well

##  The alignment process and what it produces

- we run `muscle` profile mode on a given sequence (align it to the master MSA without disturbing its columns)
- this was shown to be noisy and I will now switch to calling substitutions against the uniprot accession where it is specified.

A **substitution** is when the canonical sequence has a residue at a master alignment column but the residue differs from the column's consensus character. 
This is a genetic difference — the construct sequence differs from what the family typically has at that position.

A **deletion** is when the canonical sequence has a gap at a master alignment column that has a consensus character. This means the construct is missing a residue that most family members have.

An **insertion** is when MUSCLE introduces a new column to accommodate a residue in the canonical sequence that has no corresponding position in the master alignment. The residue exists in the canonical sequence but doesn't align to any existing family position.



##  What is actually called as a deletion — the tail problem

Here is the most important known bug and the one that is most visible in the output right now.

The disordered C-terminal tails of tubulin are frequently unresolved in crystal structures. Different constructs have different tail lengths — some are truncated by design, some are simply not resolved. When we align a canonical sequence that ends at, say, residue 430, against a master alignment whose consensus extends to position 452 (because some reference sequences have long tails), positions 431-452 in the master alignment are gaps in the canonical sequence. The current code calls each of those gap positions as a deletion variant.

- solution: mark _contiguous gap runs_ at the N/C termini (edge ~15 residues) as as TERMINAL_TRUNCATIONS rather than deletions...


##  Confidence tiers and what to do in gray areas

Not all annotations carry equal confidence, and downstream users need to know this. We propose tracking confidence through four tiers.

**Tier 1**: Reviewed UniProt accession available, strong HMM classification score, canonical sequence matches UniProt reference within a few residues. Variant calls are highly reliable. Binding site mappings are reliable.

**Tier 2**: Reviewed accession available but canonical sequence deviates noticeably from UniProt reference — possibly an engineered construct, possible sequencing discrepancy, possible genuine strain variant. Variants are called against the UniProt reference but the mismatch between PDB canonical and UniProt reference is itself flagged as requiring interpretation.

**Tier 3**: Unreviewed UniProt accession, or no accession but the entity classifies cleanly to a family and a nearest-neighbor isotype can be assigned by sequence identity. Variants are called against the nearest neighbor but flagged as lower confidence. The isotype assignment itself is uncertain.

**Tier 4**: Weak HMM classification, or the entity classifies to a family but the canonical sequence is far from any known isotype. Delta/epsilon currently mostly fall here. Annotation proceeds but with explicit warnings that the reference framework is unreliable.

For the specific case of non-vertebrate organisms: a yeast alpha tubulin or a Toxoplasma beta tubulin will often map cleanly to the master alignment in the conserved core and produce sensible substitution calls relative to the vertebrate consensus. But those substitutions are not mutations — they are phylogenetic differences expected from their evolutionary distance. Whether to include non-vertebrate sequences in the annotation output, and whether to label their "substitutions" differently, is a policy decision that depends on the use case.

---

##  What Morrissette's work gives us and where we differ

Their universal tubulin numbering system is directly comparable to our master alignment coordinate scheme. They built it by aligning 90 alpha, 90 beta, and 21 gamma sequences, cleaning out uncommon insertion columns, and numbering the resulting consensus. Their database entries — mutations, PTMs, interaction sites — are all indexed to UTN positions.

We can build a static translation table between UTN and our master alignment positions by profile-aligning their consensus sequence into our MSA. This is already partly implemented. Once that table exists, every entry in their database is directly importable into our coordinate system.

Their coverage is complementary to ours in specific ways: they have extensive curation of literature mutations, PTMs from mass spectrometry databases, and 6Å interaction environments from representative structures. We have every PDB structure computationally processed, per-chain coordinate mappings, and ligand binding sites from proximity data. Neither alone is complete.

One genuine advantage of our system: they truncate at 440 residues and explicitly exclude the CTT. If we can get the tail noise problem under control, we can index PTMs and variants in the tail — glutamylation and glycylation sites, detyrosination — which their database doesn't cover at all.

---

## Questions for Carsten


- Where are these sequences coming from that were shared a few months ago?

- Are there any positions in alpha or beta tubulin where the "wild-type" residue is genuinely different between isotypes in a way that is functionally documented — positions that define isotype identity mechanistically, not just statistically? We want to make sure those positions are correctly handled in our isotype-aware variant calling rather than flagged as mutations.

- Within vertebrates, how much genuine sequence variation exists within a single isotype across species? If we align human TUBA1A against mouse TUBA1A they're essentially identical, but where does that stop? At what phylogenetic distance does a "same isotype" comparison start to become unreliable?

- For gamma tubulin: our current MSA appears to contain non-vertebrate sequences that introduce 19 columns where all vertebrate sequences have gaps. Are these columns biologically meaningful — do they correspond to known functional regions in non-vertebrate gamma, like the specialized contacts in the gamma-TuRC lock washer — or are they just evolutionary noise? This determines whether we should mask them or keep them.

- Internal deletions in vertebrates: Are there any known genuine internal deletions in vertebrate tubulins — not tail truncations, but actual internal gaps relative to the family consensus — that are biologically meaningful? Or should we treat any internal deletion call in a vertebrate chain as a candidate artifact first?

- Gamma tub insertions: The Morrissette paper mentions that gamma tubulin has conserved insertions and deletions relative to alpha and beta that may contribute to specialized interactions. Are these insertions large enough that they would show up in our system as blocks of insertion variants in every gamma structure we process, or are they already absorbed into the gamma-specific MSA?

- Binding site definition: For binding site annotations: we define a binding site as all residues within a distance cutoff of a ligand in a structure. Is that the right conceptual unit, or do you think more in terms of defined pocket regions that should be annotated as units regardless of whether a particular structure has a ligand bound there?