# input tf -> tf count matrix -> tf energy matrix -> binding locations (genes) with affinities -> output csv (tf, tg, aff)
import json
import math
import subprocess
import pandas as pd
import random

def create_tfbs_file(gene_arr, promo_file, pcsm_database_file):
    promo_len = 100
    bind_threshold = '0.3'
    shift_len = '1'
    # lam =
    gc = '0.5'

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

    def promo_loc_to_gene(loc,key):
        index = math.floor(int(loc) / promo_len)
        return key[index]

    def run_ciiider(sequence_file, pwm_file, tf_len):
        
        config_o = f'[General]\nSTARTPOINT=1\nENDPOINT=1\n\n[Scan]\nGENELISTFILENAME={open('')}\nREFERENCEFASTA=/workspaces/PROKODE-DOCKER/src/inputs/genome.txt\nGENELOOKUPMANAGER={promo_file}\nMATRIXFILE={pwm_file}\nGENESCANRESULTS=tfsb_out.txt'
        open('config.ini', 'w').write(config_o)
        # java -jar /CiiiDER/CiiiDER_TFMs/CiiiDER.jar –n config.ini
        subprocess.run(['java','-jar','/CiiiDER/CiiiDER_TFMs/CiiiDER.jar', '-n', 'config.ini'])
        
        output_file = 'trap_results.gff'
        # annotate -s promoters.fa -psem psem.psem -t 0.1 -g 0.5 -R 21 -S 1 -o trap_results.txt 
        subprocess.run(['annotate', '-s', sequence_file, '--psem', pwm_file, '-t', bind_threshold, '-g', gc, '-R', tf_len, '-S', shift_len, '-o', output_file])
        return output_file
        
    def query_pcsm_data(gene):
        gene_found = False
        arr = ''
        tf_len = 0
        copies_found = 0
        with open(pcsm_database_file, 'r') as f:
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

    key = json.load(open('promo_key.json', 'r'))

    results = 'tf,tg,kd\n'
    for gene in gene_arr:
        print(gene)
        pscm, tf_len = query_pcsm_data(gene) #file path
        if pscm == -1:
            continue
        tf = gene

        trap_file = run_ciiider(promo_file, pscm, tf_len)
        trap_df = open(trap_file, 'r').readlines()
        trap_df = [line.replace('\n','').split('\t') for line in trap_df]
        for row in trap_df[2:]:
            tg = promo_loc_to_gene(row[3], key)
            aff =  1.0/float(row[5])
            res = f'{tf},{tg},{round(aff,5)}\n'
            print(res)
            results += res
    results = clean(results)
    open('tfbs.csv', 'w').write(results)

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

p_df = pd.read_csv('../inputs/annotation.csv')
gene_arr = p_df['geneid'].to_list()
create_decay_file(gene_arr)
create_tfbs_file(gene_arr,'promoters.fa','pcsm_data.meme')
add_betas('tfbs.csv')