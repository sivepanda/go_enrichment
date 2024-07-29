import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import requests
from io import StringIO
from Bio import Entrez

Entrez.email = "your.email@example.com"  # Replace with your email

class GOEnrichment:
    def __init__(self):
        self.go_annotations = self.fetch_go_annotations()
        self.gene_chromosome_map = self.fetch_gene_chromosome_map()
        self.data = pd.merge(self.go_annotations, self.gene_chromosome_map, on='GeneID')
    
    def fetch_go_annotations(self):
        # Fetch GO annotations from a public source
        url = 'https://geneontology.org/docs/download-go-annotations/'
        response = requests.get(url)
        response.raise_for_status()
        go_data = StringIO(response.text)
        go_annotations = pd.read_csv(go_data, sep='\t', comment='!', header=None, names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'Synonym', 'DB_Object_Type', 'Taxon_ID', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID'])
        go_annotations = go_annotations[['DB_Object_ID', 'GO_ID']]
        go_annotations.columns = ['GeneID', 'GO_ID']
        return go_annotations

    def fetch_gene_chromosome_map(self):
        # Fetch chromosome info using Ensembl REST API
        gene_ids = self.go_annotations['GeneID'].unique()
        base_url = "https://rest.ensembl.org"

        gene_chromosome_map = []

        for gene_id in gene_ids:
            ext = f"/lookup/id/{gene_id}?content-type=application/json"
            r = requests.get(base_url + ext, headers={"Content-Type": "application/json"})

            if not r.ok:
                continue

            data = r.json()
            if 'seq_region_name' in data:
                gene_chromosome_map.append({'GeneID': gene_id, 'Chromosome': data['seq_region_name']})

        return pd.DataFrame(gene_chromosome_map)
    
    def perform_enrichment(self):
        results = []

        go_ids = self.data['GO_ID'].unique()
        chromosomes = self.data['Chromosome'].unique()

        for go_id in go_ids:
            go_genes = self.data[self.data['GO_ID'] == go_id]['GeneID'].unique()
            
            for chromosome in chromosomes:
                chrom_genes = self.data[self.data['Chromosome'] == chromosome]['GeneID'].unique()
                
                overlap = len(set(go_genes) & set(chrom_genes))
                go_only = len(go_genes) - overlap
                chrom_only = len(chrom_genes) - overlap
                neither = len(self.data['GeneID'].unique()) - (overlap + go_only + chrom_only)

                contingency_table = [[overlap, go_only], [chrom_only, neither]]
                _, p_value = fisher_exact(contingency_table)
                
                results.append({
                    'GO_ID': go_id,
                    'Chromosome': chromosome,
                    'Overlap': overlap,
                    'GO_Only': go_only,
                    'Chromosome_Only': chrom_only,
                    'Neither': neither,
                    'P_Value': p_value
                })

        results_df = pd.DataFrame(results)

        # Adjust p-values for multiple testing
        results_df['Adjusted_P_Value'] = multipletests(results_df['P_Value'], method='fdr_bh')[1]

        return results_df

def main():
    enrichment = GOEnrichment()
    results = enrichment.perform_enrichment()
    results.to_csv('enrichment_results.csv', index=False)
    print("Enrichment analysis completed. Results saved to 'enrichment_results.csv'.")

if __name__ == "__main__":
    main()
