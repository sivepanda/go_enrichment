import pandas as pd
from go_enrichment.get_go import get_go_data
from scipy.stats import fisher_exact
from go_enrichment.utils import generate_contingency_table
from tqdm import tqdm

def enrichment_analysis(go_data):
    max = 10 
    go_progress_bar = tqdm(total=max, desc=f'Getting GO Data')
    
    results = []
    go_dat = []
    go_categories = set()
    chromosomes = set()
    n = 0
    
    # Collect all unique GO categories and chromosomes from input data stream
    for hit in go_data:
        go_progress_bar.update(1)
        n += 1
        if n == max: 
            go_data.close()
            break
        
        go_dat.append(hit)
        
        if 'genomic_pos_hg19' in hit:
            if len(hit['genomic_pos_hg19']) > 1:
                chromosomes.add(hit['genomic_pos_hg19'][0]['chr'])
            else:
                chromosomes.add(hit['genomic_pos_hg19']['chr'])


        for aspect in hit['go']:
            if isinstance(hit['go'][aspect], list):
                for item in hit['go'][aspect]:
                    go_categories.add(item['id'])
            else:
                go_categories.add(hit['go'][aspect]['id'])


    print("Generating Contingency tables")
    print(len(go_categories))
    for category in go_categories:
        progress_bar = tqdm(total=len(chromosomes), desc=f'Generating Contingency Tables for GO Category {category}')
        for chr in chromosomes:
            table = generate_contingency_table(go_dat, category, chr)
            progress_bar.update(1)
            oddsratio, p_value = fisher_exact(table)
            results.append({
                'go_category': category,
                'chromosome': chr,
                'odds_ratio': oddsratio,
                'p_value': p_value
            })
    
    return pd.DataFrame(results)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Perform GO enrichment analysis on chromosome data.')
    parser.add_argument('output_file', help='Path to save the enrichment analysis results CSV file.')

    args = parser.parse_args()

    go_data = get_go_data()
    results = enrichment_analysis(go_data)
    results.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()

