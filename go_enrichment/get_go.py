import mygene

# Initialize MyGene.info
mg = mygene.MyGeneInfo()

def get_go_data():
    query = 'go:*'
    fields = 'go,genomic_pos_hg19.chr'
    all_genes = mg.query(query, fields=fields, fetch_all=True)

    return all_genes
