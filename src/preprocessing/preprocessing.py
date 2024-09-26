import numpy as np
import pandas as pd
from Bio import motifs
import subprocess
import csv

temperature = 298 # kelvin

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

def get_tfbs_rowarr(line, add_betas):
    # line headers: ['tg','_','tf','tf_matrixid','start','end','strand','prescore','score','seq']
    tfbs_row = [line[2], line[0]]
    score = float(line[8])
    tfbs_row.append(score)

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
            tfbs_file.write('tf,tg,score,beta\n')
        else:
            tfbs_file.write('tf,tg,score\n')
        writer_object = csv.writer(tfbs_file)

        with open(c_output_loc, "r") as csvfile:
            prev = ''
            prev_fullline = ''
            prev_activated = False
            for row in csvfile:
                if len(row)<12:
                    row = row.replace("\n", '')
                    prev = row
                    prev_activated = True
                else:
                    if prev_activated:
                        line = prev + row
                    else:
                        line = row
                    
                    line_arr = line.split(',')
                    a = line_arr[2]
                    b = prev_fullline[2] if type(prev_fullline)==list else ''

                    if a == b:
                        prev_fullline = line_arr
                        continue
                    else: 
                        tfbs_row = get_tfbs_rowarr(line_arr, add_betas)
                        writer_object.writerow(tfbs_row)

                        prev_fullline = line_arr
                        continue
                
        return tfbs_loc

def create_promoterf(prokode_dir, genome_loc, operon_loc):
    # config files
    genome = open(genome_loc, 'r').read()
    genome = genome[genome.find('\n'):].replace('\n','')
    operons_df = pd.read_csv(operon_loc, sep='\t')
    print(operons_df)

    # write promoter region surrounding each gene's start loc to promoters.fa
    promoters_loc = prokode_dir + '/src/preprocessing/promoters.fa'    
    with open(promoters_loc, 'a') as promoters_file:
        for i in operons_df.index:
            operon = operons_df.loc[i, 'operonid']
            start = int(operons_df.loc[i, 'start'])
            promo_seq = genome[start-150:start+50]
            promoters_file.write(f'>{operon}\n{promo_seq}\n')

    return promoters_loc

def preprocessing_main(prokode_dir, genome_loc, annotation_loc, operons_loc, pfm_database_loc, CiiiDER_jar_loc, CiiiDER_thresh, add_betas=False):
    # create promoters.fa
    print('\t\t...creating promoter file')
    promoters_loc = create_promoterf(prokode_dir, genome_loc, operons_loc)

    # config motif matrix file for only tfs within genome
    print('\t\t...config tf motif file')
    motif_matrix_loc = config_tfmotifs(prokode_dir, pfm_database_loc, annotation_loc)

    # configer into tf binding site csv
    tfbs_loc = create_tfbs_through_CiiiDER(prokode_dir, CiiiDER_jar_loc, promoters_loc,motif_matrix_loc, CiiiDER_thresh, add_betas)

    return tfbs_loc
