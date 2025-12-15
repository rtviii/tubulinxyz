#!/usr/bin/env python3
"""
Example: classify a tubulin chain from a PDB structure.
"""
from pprint import pprint
from lib.hmm.tubulin_classifier import TubulinClassifier

# Initialize once (loads all HMMs)
classifier = TubulinClassifier(bitscore_threshold=50.0)

# Print HMM info
print("Loaded HMMs:")
for family, info in classifier.hmm_info().items():
    print(f"  {family}: M={info['M']}, nseq={info['nseq']}")

# Example sequence (alpha tubulin from 7SJ7)
test_sequence = """
MRECISIHVGQAGVQIGNACWELYCLEHGIQPDGQMPSDKTIGGGDDSFNTFFSETGAGKHVPR
AVFVDLEPTVIDEVRTGTYRQLFHPEQLITGKEDAANNYARGHYTIGKEIIDLVLDRIRKLADQ
CTGLQGFLVFHSFGGGTGSGFTSLLMERLSVDYGKKSKLEFAIYPAPQVSTAVVEPYNSILTTH
TTLEHSDCAFMVDNEAIYDICRRNLDIERPTYTNLNRLIGQIVSSITASLRFDGALNVDLTEFQ
TNLVPYPRIHFPLATYAPVISAEKAYHEQLSVAEITNACFEPANQMVKCDPRHGKYMACCLLY
"""

# Classify with verbose output
result = classifier.classify_verbose(test_sequence, "7SJ7", "A")

# Access programmatically
print(f"\nAssigned family: {result.assigned_family}")
print(f"Best score: {result.best_hit.score if result.best_hit else 'N/A'}")
print(f"\nAs dict:")
pprint(result.to_dict())