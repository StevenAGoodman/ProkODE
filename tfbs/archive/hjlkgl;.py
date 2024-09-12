# input tf -> tf count matrix -> tf energy matrix -> binding locations (genes) with affinities -> output csv (tf, tg, aff)
import json
import math
import subprocess
import pandas as pd
import random

def config_promo():
    gene_start_df = pd.read_csv('/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv')
    with open('/workspaces/PROKODE-DOCKER/src/inputs/genome.txt', 'r') as gf:
        genome = gf.read() ## LEARN HOW TO READ GENOME MORE EFFICIENTLY (KW: STRING SEEK)

    key = []
    poutput = ''
    goutput = ''
    for i in gene_start_df.index:
        gene = gene_start_df.loc[i, 'geneid']
        start = gene_start_df.loc[i, 'start']
        promo_seq = genome[start-150:start+50]
        poutput += f'>{gene}\n{promo_seq}\n'

    promoter_file = '/workspaces/PROKODE-DOCKER/src/preprocessing/promoters.fa'
    with open(promoter_file, 'w') as file:
        file.write(poutput)
    return promoter_file
        
def create_tfbs_file(gene_list_file, promo_file, pscm_database_file):
    gene_arr = open(gene_list_file, 'r').read().split('\n')

    def clean(csv_str):
        return_str = ""
        prev = ''
        for row in csv_str.split('\n'):
            a = row.find(',')
            b = prev.find(',')
            c = row[:(a+1 + row[a+1:].find(','))]
            d = prev[:(b+1 + prev[b+1:].find(','))]

            if c == d:
                aff = sum([float(row[(a+2 + row[a+1:].find(',')):]), float(prev[(b+2 + prev[b+1:].find(',')):])])/2.0
                prev = f'{c},{aff}'
            else:
                return_str += prev + '\n'
                prev = row
        return return_str

    def run_ciiider(gene_list_file, pfm_file,):
        output_file = '/workspaces/PROKODE-DOCKER/src/preprocessing/ciiider_results.txt'
        config_o = f'[General]\nSTARTPOINT=1\nENDPOINT=1\n\n[Scan]\nGENELISTFILENAME={gene_list_file}\nREFERENCEFASTA=/workspaces/PROKODE-DOCKER/src/inputs/genome.fa\nGENELOOKUPMANAGER={promo_file}\nMATRIXFILE={pwm_file}\nGENESCANRESULTS={}'
        open('/workspaces/PROKODE-DOCKER/src/preprocessing/config.ini', 'w').write(config_o)
        # java -jar /CiiiDER/CiiiDER_TFMs/CiiiDER.jar –n config.ini
        subprocess.run(['java','-jar','/CiiiDER/CiiiDER_TFMs/CiiiDER.jar', '-n', 'config.ini'])
        
        return output_file
        
    def query_pcsm_data(gene):
        gene_found = False
        arr = ''
        tf_len = 0
        copies_found = 0
        with open(pscm_database_file, 'r') as f:
            for _, line in enumerate(f):
                if line.casefold().find(gene.casefold()) >= 0:
                    gene_found = True
                    copies_found+=1
                if gene_found == True:
                    if line == '\n':
                        gene_found = False
                        arr += line
                    else:
                        arr += line
                        tf_len += 1

        fname = 'pscm.pscm'
        if arr == '':
            fname = -1
        else:
            with open(fname, 'w') as f:
                f.write(arr)
            tf_len = int(tf_len/copies_found)
        return fname, f'{tf_len}'

    promo_file = config_promo()

    # ciiider_results_f = run_ciiider(promo_file, pscm_database_file)

    # results = 'tf,tg,kd\n'
    # for gene in gene_arr:
    #     print(gene)
    #     pscm, tf_len = query_pcsm_data(gene) #file path
    #     if pscm == -1:
    #         continue
    #     tf = gene

    #     ciiider_file = run_ciiider(gene_list_file, pscm, tf_len)
    #     trap_df = open(ciiider_file, 'r').readlines()
    #     trap_df = [line.replace('\n','').split('\t') for line in trap_df]
    #     for row in trap_df[2:]:
    #         tg = promo_loc_to_gene(row[3], key)
    #         aff =  1.0/float(row[5])
    #         res = f'{tf},{tg},{round(aff,5)}\n'
    #         print(res)
    #         results += res
    # results = clean(results)
    # open('tfbs.csv', 'w').write(results)

create_tfbs_file('gene_list.fa','promoters.fa','pfmdb.txt')

def create_decay_file(gene_arr):
    output = 'gene,decay\n'
    for gene in gene_arr:
        # gene_seq = 
        decay = random.uniform(0,0.001) # DegScore(gene_seq, ('.' * len(gene_seq))).est_k_deg
        output += f'{gene},{round(decay,5)}\n'
    
    open('../decay_rates.csv', 'w').write(output)

def add_betas(csv):
    csv_df = pd.read_csv(csv)

    col = []
    for tf in csv_df['tf'].values:
        seq = tf #jal;dfj
        beta = random.uniform(0.0,100.0)
        col.append(round(beta,4))

    csv_df.insert(3,'beta', col,allow_duplicates=True)
    csv_df.to_csv(open('../tfbs.csv', 'w'))

