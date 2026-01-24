
# Quick test script - save as test_classify.py and run it
from lib.etl.molstar_bridge import run_molstar_extraction
from lib.etl.classification import classify_sequence, get_classifier
from lib.types import ObservedSequenceData, ObservedResidue
from pathlib import Path
from api.config import PROJECT_ROOT

# Test with a known tubulin structure
rcsb_id = "5JCO"  # A well-known tubulin structure

# First, let's just test the classifier directly with a known alpha-tubulin sequence fragment
test_seq = "MRECISIHVGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGKHVPRAVFVDLEPTVIDEVRTGTYRQLFHPEQLITGKEDAANNYARGHYTIGKEIIDLVLDRIRKLADQCTGLQGFLVFHSFGGGTGSGFTSLLMERLSVDYGKKSKLEFSIYPAPQVSTAVVEPYNSILTTHTTLEHSDCAFMVDNEAIYDICRRNLDIERPTYTNLNRLIGQIVSSITASLRFDGALNVDLTEFQTNLVPYPRIHFPLATYAPVISAEKAYHEQLSVAEITNACFEPANQMVKCDPRHGKYMACCLLYRGDVVPKDVNAAIATIKTKRTIQFVDWCPTGFKVGINYQPPTVVPGGDLAKVQRAVCMLSNTTAIAEAWARLDHKFDLMYAKRAFVHWYVGEGMEEGEFSEAREDMAALEKDYEEVGIDSYEDEDEGEE"

observed = ObservedSequenceData(
    auth_asym_id="A",
    entity_id="1",
    residues=[
        ObservedResidue(auth_seq_id=i+1, label_seq_id=i+1, comp_id="ALA", one_letter=aa)
        for i, aa in enumerate(test_seq)
    ]
)

print(f"Testing classification with sequence of length {len(test_seq)}")
print(f"First 50 chars: {test_seq[:50]}")

# Get classifier and run
classifier = get_classifier()
print(f"Classifier has {len(classifier.hmms)} HMMs loaded")
print(f"HMM families: {[f.value for f in classifier.hmms.keys()]}")

# Run classification
result = classifier.classify(test_seq, "TEST", "A")
print(f"\nClassification result:")
print(f"  Assigned: {result.assigned_family}")
print(f"  Best hit: {result.best_hit}")
print(f"  All hits ({len(result.hits)}):")
for hit in sorted(result.hits, key=lambda h: h.score, reverse=True)[:5]:
    print(f"    {hit.family.value}: score={hit.score:.1f}, evalue={hit.evalue:.2e}")
