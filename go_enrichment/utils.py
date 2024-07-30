import json

def generate_contingency_table(go_data, go_category, chromosome):
    chr_with_category = 0
    chr_without_category = 0
    non_chr_with_category = 0
    non_chr_without_category = 0
    x = 0
    for hit in go_data:
        if 'genomic_pos_hg19' in hit:
            if len(hit['genomic_pos_hg19']) > 1:
                chr = (hit['genomic_pos_hg19'][0]['chr'])
            else:
                chr = (hit['genomic_pos_hg19']['chr'])
        
        has_category = any(go['id'] == go_category for aspect in hit['go'] for go in hit['go'][aspect] if isinstance(hit['go'][aspect], list))
        
        if chr == chromosome:
            if has_category:
                chr_with_category += 1
            else:
                chr_without_category += 1
        else:
            if has_category:
                non_chr_with_category += 1
            else:
                non_chr_without_category += 1
    
    return [
        [chr_with_category, chr_without_category],
        [non_chr_with_category, non_chr_without_category]
    ]

