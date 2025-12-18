from Bio import SeqIO

def analyze_cluster_composition(
    raw_fasta: str,
    clstr_file: str,
    top_n_clusters: int = 5,
):
    """
    Look inside the largest clusters to see what we're collapsing.
    Are we losing meaningful isotype diversity or just redundant orthologs?
    """
    from collections import defaultdict
    import re
    
    # Load all sequences with metadata
    records = {r.id.split("|")[1]: r for r in SeqIO.parse(raw_fasta, "fasta")}
    
    # Parse clusters
    clusters = []
    current_cluster = {"id": None, "members": []}
    
    with open(clstr_file) as f:
        for line in f:
            if line.startswith(">Cluster"):
                if current_cluster["members"]:
                    clusters.append(current_cluster)
                current_cluster = {"id": line.strip(), "members": []}
            else:
                # Extract accession
                match = re.search(r'>([^|]+\|)?([A-Z0-9]+)', line)
                if match:
                    acc = match.group(2)
                    is_rep = line.strip().endswith("*")
                    current_cluster["members"].append({"acc": acc, "is_rep": is_rep})
        
        if current_cluster["members"]:
            clusters.append(current_cluster)
    
    # Sort by size
    clusters.sort(key=lambda c: len(c["members"]), reverse=True)
    
    print(f"\nAnalyzing top {top_n_clusters} clusters:\n")
    
    for cluster in clusters[:top_n_clusters]:
        print(f"{cluster['id']} - {len(cluster['members'])} members")
        
        # Group by organism and gene name
        by_organism = defaultdict(list)
        gene_names = set()
        
        for member in cluster["members"]:
            acc = member["acc"]
            if acc in records:
                rec = records[acc]
                desc = rec.description
                
                # Extract organism
                org_match = re.search(r'OS=([^=]+)\s+OX=', desc)
                organism = org_match.group(1).strip() if org_match else "Unknown"
                
                # Extract gene name
                gene_match = re.search(r'GN=(\S+)', desc)
                gene = gene_match.group(1) if gene_match else "N/A"
                
                by_organism[organism].append(gene)
                gene_names.add(gene)
                
                if member["is_rep"]:
                    print(f"  REP: {acc} | {gene} | {organism}")
        
        # Summary
        print(f"  Organisms: {len(by_organism)}")
        print(f"  Unique gene names: {gene_names}")
        
        # Show organism breakdown
        top_orgs = sorted(by_organism.items(), key=lambda x: -len(x[1]))[:5]
        for org, genes in top_orgs:
            print(f"    {org}: {len(genes)} seqs, genes: {set(genes)}")
        
        print()


# Run on alpha
analyze_cluster_composition(
    "./raw_fastas/tubulin_alpha_raw.fasta",
    "./clustered/tubulin_alpha_nr.fasta.clstr",
    top_n_clusters=5
)