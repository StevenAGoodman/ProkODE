import numpy as np
import pandas as pd
from Bio import motifs
import subprocess
import csv

def config_tfmotifs(prokode_dir, pfm_database_loc, annotation_loc):
    motif_arr = []

    annotation_df = pd.read_csv(annotation_loc, delimiter='\t')
    gene_vector =  annotation_df['geneid'].tolist()
    record = motifs.parse(open(pfm_database_loc, 'r'), 'JASPAR')
    
    for m in record:
        if m.name in gene_vector:
            motif_arr.append(m)
        else:
            continue
    
    out_loc = prokode_dir + '/src/preprocessing/motif_matrices.txt'
    out = motifs.write(motif_arr, "JASPAR")
    open(out_loc, 'w').write(out)

    return out_loc

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

def write_to_tfbs(line, add_betas):
    # line headers: ['tg','_','tf','tf_matrixid','start','end','strand','prescore','score','seq']
    tfbs_row = [line[2], line[0]]
    kd = score_to_kd(line[8])
    tfbs_row.append(kd)

    # add betas if desired
    if add_betas:
        beta = np.random.randint(0,100)
        tfbs_row.append(beta)
    
    return tfbs_row


def create_tfbs_through_CiiiDER(prokode_dir, jar_loc, promoters_loc, matrices_loc, deficit_val, add_betas):
    # CiiiDER (https://ciiider.erc.monash.edu/) is a software that searches the promoter DNA regions of each gene with the binding motifs of each transcription factor to determine their binding sites. 
    c_output_loc = prokode_dir + '/src/preprocessing/CiiiDER_results.txt'
    
    # config config.ini
    config_o = f"""[General]
STARTPOINT = 1
ENDPOINT = 1\n
[Scan]
GENELISTFILENAME = {promoters_loc}
MATRIXFILE = {matrices_loc}
GENESCANRESULTS = {c_output_loc}
DEFICIT = {deficit_val}"""
    
    open(prokode_dir + '/src/preprocessing/config.ini', 'w').write(config_o)

    subprocess.run(['java','-jar', jar_loc, '-n', prokode_dir + '/src/preprocessing/config.ini'])

    c_output_loc = c_output_loc[:-3] + 'csv'

    # cleanup results and write to tfbs.csv
    print('\t\t...creating tf binding site file')
    tfbs_loc = prokode_dir + '/src/tfbs.csv'
    with open(tfbs_loc, 'a') as tfbs_file:
        if add_betas:
            tfbs_file.write('tf,tg,kd,beta\n')
        else:
            tfbs_file.write('tf,tg,kd\n')
        writer_object = csv.writer(tfbs_file)

        with open(c_output_loc, "r") as csvfile:
            prev = ''
            prev_fullline = ''
            prev_act = False
            for row in csvfile:
                if len(row)<10:
                    row = row.replace("\n", '')
                    prev = row
                    prev_act = True
                else:
                    if prev_act:
                        line = prev + row
                    else:
                        line = row
                    
                    line = line.split(',')
                    a = line[1]
                    b = prev_fullline[1]

                    if a == b:
                        continue
                    else: 
                        tfbs_row = write_to_tfbs(line, add_betas)
                        writer_object.writerow(tfbs_row)

                        prev_fullline = line
                
        return tfbs_loc

def create_promoterf(prokode_dir, genome_loc, annotation_loc):
    # config files
    genome = open(genome_loc, 'r').read()
    annotation_df = pd.read_csv(annotation_loc, delimiter='\t')

    # write promoter region surrounding each gene's start loc to promoters.fa
    promoter_contents = ''
    for i in annotation_df.index:
        gene = annotation_df.loc[i, 'geneid']
        start = int(annotation_df.loc[i, 'start'])
        promo_seq = genome[start-150:start+50]
        promoter_contents += f'>{gene}\n{promo_seq}\n'

    # write to promoters.fa
    promoters_loc = prokode_dir + '/src/preprocessing/promoters.fa'
    open(promoters_loc, 'w').write(promoter_contents)
    return promoters_loc

def score_to_kd(score):
    return np.exp(-1 * score)

def create_tfbs(prokode_dir, ci_results_loc):
    # import ciiider results
    ciiider_df = pd.read_csv(ci_results_loc, names=['tg','_','tf','tf_matrixid','start','end','strand','prescore','score','seq'])

    tfbs_df = ciiider_df[['tf','tg']].copy()

    kd_vals = []
    for score in ciiider_df['score'].to_list():
        kd = score_to_kd(score)
        kd_vals.append(kd)

    tfbs_df.insert(2,'Kd', kd_vals)
    tfbs_df.drop(axis=1, index=0)

    # write to tfbs.csv
    tfbs_loc = prokode_dir + '/src/tfbs.csv'
    out = tfbs_df.to_csv()
    out = clean(out)
    open(tfbs_loc, 'w').write(out)
    
    return tfbs_loc

def get_decay_rates(gene):
    return np.log(2)/300

def decay_rates_main(prokode_dir, annotation_loc):
    # read number of genes
    annotation_df = pd.read_csv(annotation_loc, delimiter='\t')

    out = 'gene,decay_rate\n'
    for gene in annotation_df['geneid'].to_list():
        decay_rate = get_decay_rates(gene)
        out += f'{gene},{decay_rate}\n'

    decay_rates_loc = prokode_dir + '/src/decay_rates.csv'
    open(decay_rates_loc, 'w').write(out)

    return decay_rates_loc

def preprocessing_main(prokode_dir, genome_loc, annotation_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, add_betas=False):
    # create promoters.fa
    print('\t\t...creating promoter file')
    promoters_loc = create_promoterf(prokode_dir, genome_loc, annotation_loc)

    # config motif matrix file for only tfs within genome
    print('\t\t...config tf motif file')
    motif_matrix_loc = config_tfmotifs(prokode_dir, pfm_database_loc, annotation_loc)

    # configer into tf binding site csv
    tfbs_loc = create_tfbs_through_CiiiDER(prokode_dir, CiiiDER_jar_loc, promoters_loc,motif_matrix_loc, CiiiDER_thresh, add_betas)

    # create decay_rates.csv
    print('\t\t...creating decay rates file')
    decay_rates_loc = decay_rates_main(prokode_dir, annotation_loc)
    
    return tfbs_loc, decay_rates_loc

# # temp start code
# preprocessing_main('/workspaces/PROKODE-DOCKER','/workspaces/PROKODE-DOCKER/src/inputs/genome.fasta','/workspaces/PROKODE-DOCKER/src/inputs/annotation.csv','/workspaces/PROKODE-DOCKER/src/preprocessing/pfmdb.txt')
