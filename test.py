from mygene import MyGeneInfo
import pyensembl

# Initialize MyGeneInfo and Ensembl data
mg = MyGeneInfo()
ensembl = pyensembl.EnsemblRelease(release=95)

# Fetch all genes on a specific chromosome, e.g., chromosome 1
genes_chr1 = ensembl.gene_ids(contig='1')

# for gene in genes_chr1:
    # print(gene)

# Get gene symbols
# gene_symbols = [gene.gene_name for gene in genes_chr1]

# Query MyGeneInfo for GO terms
gene_info = mg.querymany(genes_chr1, scopes='symbol', fields='go', species='human')

# Print GO terms for each gene
for gene in gene_info:
    print(f"Gene: {gene['query']}")
    if 'go' in gene:
        for go_category in gene['go']:
            print(f"  {go_category['term']}: {go_category['id']}")
    else:
        print("  No GO terms found")

