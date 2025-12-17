- Tubulin Classification (via hmms)
- ligand registry (figure out if possible) 

# MIP/MAP classification

# The landing page:

After everything else is done:
- purported ligand bindign sites
- mutation hotspots
- annotations 

# The structure page

- annotations registry

# The MSA page

- individual chain structural superposition

--------------------------------------------

<!-- - db lib reader implement -->
- create master alignment models and nodes
- implement the mapping on the links
- makes sure the ligands are stored as unique nodes, no duplication in the profile or DB
- implement entity_id->auth_asym_id backlinking system
- write and test in/del/mut logic  
- PTM logic
- Morisette dataset of PTMs

--------------------------------------------

1. HMM classification models: alpha, beta, main MIPs/MAPs
- for mip/maps -- record references

- add classification (major families & mip/maps)
- search non-canonical sequences in pdb for modifications

2. Standardize mutations, modifications (between morisette and structural)
- make sure the morisette stuff corresponds to the master alignment

3. Ligands interfaces (Nonpolymer)->PolymerInstance
- how to parse and how to store? 

### fend:
- msa-viewer
- interactions for the msa-viewer
- new colorschme for the molstar preset