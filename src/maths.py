import numpy as np
import pandas as pd

def eq1(beta, N_tf, Kd_tf, N_p, Kd_p, Nns):
    P_tf = N_tf / (N_tf + Kd_tf)
    E_basal = N_p / (Nns + Kd_p)
    E_PTF = beta * E_basal
    R = P_tf * E_PTF + (1 - P_tf) * E_basal
    return R

def rev_eq1(R,N_tf, Kd_tf, N_p, Kd_p, Nns):
    P_tf = N_tf / (N_tf + Kd_tf)
    E_basal = N_p / (Nns + Kd_p)
    E_PTF = (R - ((1 - P_tf) * E_basal)) / P_tf
    beta = E_PTF / E_basal
    return beta

# all the data we need for a tf-rnap-tf interaction is the prob of rnap binding at instant, N_tf, Kd_tf, N_p, Kd_p, Nns
# how does gene expression relate to prob of rnap binding at instant??? 

def subset_tf_of_genes(gene_arr):

    return tf_arr

def get_Kd(binder, bindee):
    Kd =
    return Kd

def express_to_R(exp):
    return R

def get_tgs(tf):
    return tg_arr

def main():
    df = pd.read_csv('data.csv', names=['gene', 'expression'])
    tfs = subset_tf_of_genes(df['gene'].values)

    Nns = 4600000
    p_gene = smth?
    N_p = df.loc[p_gene, 'expression']

    with open('results.csv', 'w') as file:
        file.write('tf,tg,beta')
        for tf in tfs:
            N_tf = df.loc[tf, 'expression']
            tg_arr = get_tgs(tf)
            for tg in tg_arr:
                Kd_tf = get_Kd(tf, tg)
                Kd_p = get_Kd(p_gene, tg)
                R = express_to_R()
                beta = rev_eq1(R, N_tf, Kd_tf, N_p,Kd_p,Nns)
                file.write(f'{tf},{tg},{beta}')

