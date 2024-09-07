import pandas as pd

gene_start_df = pd.read_csv('../inputs/annotation.csv')
with open('../inputs/genome.txt', 'r') as gf:
    genome = gf.read() ## LEARN HOW TO READ GENOME MORE EFFICIENTLY (KW: STRING SEEK)

key = []
poutput = ''
goutput = ''
for i in gene_start_df.index:
    gene = gene_start_df.loc[i, 'geneid']
    start = gene_start_df.loc[i, 'start']
    promo_seq = genome[start-100:start+50]
    poutput += f'>{gene}\n{promo_seq}\n'
    goutput += f'{gene}\n'

with open('promoters.fa', 'w') as file:
    file.write(poutput)
with open('gene_list.fa', 'w') as file:
    file.write(goutput)