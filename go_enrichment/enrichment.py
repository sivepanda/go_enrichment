import pandas as pd
from go_enrichment.get_go import get_go_data
from scipy.stats import fisher_exact
from go_enrichment.utils import generate_contingency_table
from tqdm import tqdm



# This program uses the mygene.info API, an example you can visit to see usual output of this is below
# https://mygene.info/v3/query?q=go:*&fields=go,genomic_pos_hg19.chr&fetch_all=TRUE

def process_chromosome(go_dat, go_info, id, chr):
    # for chr in chromosomes:
    table = generate_contingency_table(go_dat, id , chr)
    progress_bar.update(1)
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

def enrichment_analysis(go_data, num_points):
    go_progress_bar = tqdm(total=num_points, desc=f'Getting GO Data')
    
    results = []
    go_dat = []
    go_ids = set()
    chromosomes = set()
    go_info = {} # Dictionary that has the GO ID as the Key and the Term and Category within it
    n = 0
    
    # Collect all unique GO categories and chromosomes from input data stream
    for hit in go_data:
        go_progress_bar.update(1)
        n += 1
        if n == num_points: 
            go_data.close()
            break
        
        go_dat.append(hit)
        
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



    print("Generating Contingency tables")
    for id in go_ids:
        progress_bar = tqdm(total=len(chromosomes), desc=f'Generating Contingency Tables for GO Category')
         with concurrent.futures.ThreadPoolExecutor() as executor:
                for id in go_ids:
                    progress_bar = tqdm(total=len(chromosomes), desc=f'Generating Contingency Tables for GO Category {id}')
                    # Create a list of futures
                    futures = {executor.submit(process_chromosome, go_dat, go_info, id, chr): chr for chr in chromosomes}
                    for future in concurrent.futures.as_completed(futures):
                        result = future.result()
                        results.append(result)
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

