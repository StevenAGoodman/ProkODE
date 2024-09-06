# output: TF_ID,TG_ID,R,N_tf,Kd_tf,N_p,Kd_p,N_ns
import numpy as np
import pandas as pd

tfbs_file = '../src/processing/results.csv'
N_ns = 4600000
poly_gene = ''

def get_tfs(tg):
    tf_tgs = tfbs_df.loc[tfbs_df['tg'] == tg]
    tf_arr = tf_tgs['tf'].to_list()
    tf_kd_arr = tf_tgs['kd'].to_list()

    return tf_arr, tf_kd_arr # tfs that regulate tg

def get_Kd_p(tg):
    # poly or sigma?
    return Kd_p

def tf_knockout(data_files):
    
    result = ''
    for file in data_files:
        pd.read_csv()

    return result

def time_series(data_files):
    
    result = ''
    for file in data_files:
        pd.read_csv()

    return result

def chip_seq(data_files):
    
    result = ''
    for file in data_files:
        pd.read_csv()

    return result

def static(data_files):

    result = ''
    for file in data_files:
        df = pd.read_csv(file, names=['gene', 'expression']) ### FILE FORMAT ###

        total_mRNA = sum(df['expression'].values)
        N_p = df.loc(poly_gene, "expression")

        for tg_row in df.iterrows: 
            TG_ID = tg_row['gene']
            Kd_p = get_Kd_p(TG_ID)
            R = tg_row['expression'] / total_mRNA
            tf_array, tf_kd_array = get_tfs(TG_ID)

            for i in range(len(tf_array)):
                N_tf = df.loc[tf_array[i], 'expression']
                result += f'{tf_array[i]},{TG_ID},{R},{N_tf},{tf_kd_array[i]},{N_p},{Kd_p},{N_ns}\n'

    return result

tfbs_df = pd.read_csv(tfbs_file)
data_files = []
return_txt = # fucntion(data_files)
results = open('results/result.csv', 'a')
results.write(return_txt)