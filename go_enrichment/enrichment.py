import pandas as pd
from go_enrichment.get_go import get_go_data
from scipy.stats import fisher_exact
from go_enrichment.utils import generate_contingency_table
from tqdm import tqdm
import concurrent.futures
import numpy as np



# This program uses the mygene.info API, an example you can visit to see usual output of this is below
# https://mygene.info/v3/query?q=go:*&fields=go,genomic_pos_hg19.chr&fetch_all=TRUE

# Completes a Fisher Exact test on a GO Category and the chromosomes it is found in
def process_chromosome(go_dat, go_info, id, chr):
    table = generate_contingency_table(go_dat, id , chr)
    oddsratio, p_value = fisher_exact(table)
    result = ({
        'chromosome': chr,
        'go_id': id,
        'go_name': go_info[id]['term'],
        'top lvl': go_info[id]['category'],
        'GO: Y; Chr :Y': table[0][0],
        'GO: N; Chr :Y': table[0][1],
        'GO: Y; Chr :N': table[1][0],
        'GO: N; Chr :N': table[1][1],
        'odds_ratio': oddsratio,
        'p-value': p_value
    })
    return result

# Runs the enrichment analysis
def enrichment_analysis(go_data, num_points):
    possible_chr = values = [str(i) for i in range(1, 24)] + ['X', 'Y']
    go_progress_bar = tqdm(total=num_points, desc=f'Getting GO Data') # Progress bar for fetching GO data.
    
    go_info = {} # Dictionary that has the GO ID as the Key and the Term and Category within it
    go_dat = [] # To preserve the GO data after the user limit has been reached (object is destroyed when the datastream ends)
    results = []
    go_ids = set()
    chromosomes = set()
    n = 0
    
    # Collect all unique GO categories and chromosomes from input data stream
    for hit in go_data:
        go_progress_bar.update(1)
        n += 1
        
        # Exit when max is hit
        if n == num_points: 
            go_dat_np = np.array(go_dat)
            np.save('./go_data.npy', go_dat_np)
            go_data.close()
            break

        if n % 10 == 0: 
            go_dat_np = np.array(go_dat)
            np.save('./go_data.npy', go_dat_np)

        # Do not process when chromosome is not a known one 1-23, X, Y
        if 'genomic_pos_hg19' in hit:
            if len(hit['genomic_pos_hg19']) > 1:
                if not ( str(hit['genomic_pos_hg19'][0]['chr']) in possible_chr):
                    if ( str(hit['genomic_pos_hg19'][1]['chr']) in possible_chr ):
                        hit['genomic_pos_hg19'][0]['chr'] = hit['genomic_pos_hg19'][1]['chr']
                    else:
                        continue
            else:
                if not ( str(hit['genomic_pos_hg19']['chr']) in possible_chr ):
                    continue
        
        go_dat.append(hit)

    

    print(len(go_dat))
    for hit in go_dat:
        if 'genomic_pos_hg19' in hit:
            if len(hit['genomic_pos_hg19']) > 1:
                chromosomes.add(hit['genomic_pos_hg19'][0]['chr'])
            else:
                chromosomes.add(hit['genomic_pos_hg19']['chr'])

        # So that both the BP and MP are queried for each response
        for aspect in hit['go']:
            if isinstance(hit['go'][aspect], list):
                for item in hit['go'][aspect]:
                    go_ids.add(item['id'])
                    go_info[item['id']] = { 'category': aspect , 'term': item['term'] }
            else:
                go_ids.add(hit['go'][aspect]['id'])
                go_info[hit['go'][aspect]['id']] = { 'category': aspect, 'term':  hit['go'][aspect]['term'] }


    num_processed = 0
    print("Generating Contingency tables")
    with concurrent.futures.ThreadPoolExecutor() as executor:
       for id in go_ids:
           progress_bar = tqdm(total=len(chromosomes), desc=f'Generating Contingency Tables for GO Category {id}')
           # Create a list of futures
           futures = {executor.submit(process_chromosome, go_dat, go_info, id, chr): chr for chr in chromosomes}
           for future in concurrent.futures.as_completed(futures):
                result = future.result()
                results.append(result)
                if num_processed % 100 == 0:
                    results_np = np.array(results)
                    np.save('./results.npy', results_np)
                progress_bar.update(1)   
    
    return pd.DataFrame(results)

def main():
    import argparse

    parser = argparse.ArgumentParser(description='Perform GO enrichment analysis on chromosome data.')
    parser.add_argument('output_file', help='Path to save the enrichment analysis results CSV file.')
    parser.add_argument('num_points', help='The number of GO IDs you would like to query.')

    args = parser.parse_args()

    go_data = get_go_data()
    results = enrichment_analysis(go_data, int(args.num_points))
    results.to_csv(args.output_file, index=False)

if __name__ == "__main__":
    main()

