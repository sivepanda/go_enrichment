import requests
import pandas as pd
from goatools.base import get_godag
from goatools.anno.genetogo_reader import Gene2GoReader

def fetch_go_annotations():
    # Using goatools to fetch GO annotations
    godag = get_godag("http://current.geneontology.org/ontology/go-basic.obo", optional_attrs="relationship")
    g2g = Gene2GoReader("http://geneontology.org/gene-associations/goa_human.gaf.gz", taxids=[9606])
    
    gene2go = g2g.get_ns2nts(BP)  # Biological Process GO annotations
    go_annotations = []

    for geneid, go_terms in gene2go.items():
        for go_term in go_terms:
            go_annotations.append({GeneID: geneid, GO_ID: go_term.GO_id})

    return pd.DataFrame(go_annotations)

def fetch_gene_chromosome_map():
    # Fetch chromosome info using Ensembl REST API
    # Example using a small set of gene IDs for demonstration purposes
    gene_ids = ["ENSG00000139618", "ENSG00000157764", "ENSG00000166710"]
    base_url = "https://rest.ensembl.org"

    gene_chromosome_map = []

    for gene_id in gene_ids:
        ext = f"/lookup/id/{gene_id}?content-type=application/json"
        r = requests.get(base_url + ext, headers={"Content-Type": "application/json"})

        if not r.ok:
            continue

        data = r.json()
        gene_chromosome_map.append({GeneID: gene_id, Chromosome: data[seq_region_name]})

    return pd.DataFrame(gene_chromosome_map)
