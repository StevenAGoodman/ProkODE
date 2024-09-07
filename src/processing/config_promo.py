import pandas as pd
import json

gene_start_df = pd.read_csv('../inputs/annotation.csv')
with open('../inputs/genome.txt', 'r') as gf:
    genome = gf.read() ## LEARN HOW TO READ GENOME MORE EFFICIENTLY (KW: STRING SEEK)

key = []
output = '>divide position by 100 and round to get index of key (in config_promo.py)\n'
for i in gene_start_df.index:
    gene = gene_start_df.loc[i, 'geneid']
    start = gene_start_df.loc[i, 'start']
    seq = genome[start-70:start+30]
    output += f'{seq}\n'
    key.append(gene)

with open('promoters.fa', 'w') as file:
    file.write(output)
with open('promo_key.json', 'w') as file:
    data = f'{key}'.replace("'", '"')
    file.write(data)