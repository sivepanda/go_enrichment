import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from .fetch_data import fetch_go_annotations, fetch_gene_chromosome_map

class GOEnrichment:
    def __init__(self):
        self.go_annotations = fetch_go_annotations()
        self.gene_chromosome_map = fetch_gene_chromosome_map()
        self.data = pd.merge(self.go_annotations, self.gene_chromosome_map, on=GeneID)
    
    def perform_enrichment(self):
        results = []

        go_ids = self.data[GO_ID].unique()
        chromosomes = self.data[Chromosome].unique()

        for go_id in go_ids:
            go_genes = self.data[self.data[GO_ID] == go_id][GeneID].unique()
            
            for chromosome in chromosomes:
                chrom_genes = self.data[self.data[Chromosome] == chromosome][GeneID].unique()
                
                overlap = len(set(go_genes) & set(chrom_genes))
                go_only = len(go_genes) - overlap
                chrom_only = len(chrom_genes) - overlap
                neither = len(self.data[GeneID].unique()) - (overlap + go_only + chrom_only)

                contingency_table = [[overlap, go_only], [chrom_only, neither]]
                _, p_value = fisher_exact(contingency_table)
                
                results.append({
                    GO_ID: go_id,
                    Chromosome: chromosome,
                    Overlap: overlap,
                    GO_Only: go_only,
                    Chromosome_Only: chrom_only,
                    Neither: neither,
                    P_Value: p_value
                })

        results_df = pd.DataFrame(results)

        # Adjust p-values for multiple testing
        results_df[Adjusted_P_Value] = multipletests(results_df[P_Value], method=fdr_bh)[1]

        return results_df

def main():
    enrichment = GOEnrichment()
    results = enrichment.perform_enrichment()
    results.to_csv(enrichment_results.csv, index=False)
    print("Enrichment analysis completed. Results saved to enrichment_results.csv.")

