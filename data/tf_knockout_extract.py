import numpy as np
import pandas as pd

def eq1(R,N_tf, Kd_tf, N_p, Kd_p, Nns):
    P_tf = N_tf / (N_tf + Kd_tf)
    E_basal = N_p / (Nns + Kd_p)
    E_PTF = (R - ((1 - P_tf) * E_basal)) / P_tf
    beta = E_PTF / E_basal
    return beta

import subprocess
import pandas as pd
import numpy as np
from scipy.optimize import minimize

### FURTHER READING FOR SOME PROBLEMS
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5910820/

def getAlphaBeta(expression_nonzero, expression_atzero, num_tfs, num_rnap_ctrl, num_rnap_treat, num_totalgenes, tf_bindAff, rnap_bindAff):
    m = (num_tfs / num_totalgenes)
    b_ctrl = (num_rnap_ctrl / num_totalgenes)
    # b_treat = (num_rnap_treat / num_totalgenes) ## perhaps look into further?
    delta_tf = np.e ** (-1 * tf_bindAff)
    delta_p = np.e ** (-1 * rnap_bindAff)

    beta = ((((expression_atzero - 1) * (b_ctrl * delta_p + 1)) / (m * delta_tf)) - 1) / (b_ctrl * delta_p)
    alpha = (expression_nonzero - 1) / (delta_tf * beta * m)
    return alpha, beta

def get_tf_bindaff(tg_loc, tf_bindlocs):
    chr = tg_loc[0]
    #narrow tf_bindlocs by chrm
    tf_bindlocs = tf_bindlocs.loc[tf_bindlocs['chr']==chr]
    tf_bindlocs = tf_bindlocs[tf_bindlocs['start']==min(tf_bindlocs['start'].tolist(), key=lambda x:abs(x-tg_loc[1]))]
    tf_bindscore = tf_bindlocs.iloc[0]['score']
    return tf_bindscore

def tg_main(tg_arr, tf_bindlocs_bindscores, tf_amnt, rnap_amnt_ctrl, rnap_amnt_treat, rnap_bindaff):
    target_gene = tg_arr[0]
    ctrl_expression = tg_arr[2]
    treat_expression = tg_arr[3]

    tg_loc = hg19_promoter_df.loc[hg19_promoter_df['name'].str.contains(f'{target_gene}_.')]
    if tg_loc.empty:
        return target_gene, 'NaN','NaN'
    tg_loc = tg_loc.iloc[0].values.flatten().tolist()
    tf_bindaff = get_tf_bindaff(tg_loc, tf_bindlocs_bindscores)  * 0.001

    err_func = lambda x: (30.3379 - ((1 + x[0] * x[1] * (100 / 3.1*10**9) * np.e ** -0.815) / (1 + (1 + x[1] * (100 * np.e ** -0.7)) * np.e ** -0.815 / 3.1*10**9)))**2
    getAlphaBeta = minimize(err_func, [10,123], tol=1e-6)
    print(getAlphaBeta)
    # alpha, beta = getAlphaBeta(ctrl_expression,treat_expression,tf_amnt,rnap_amnt_ctrl,rnap_amnt_treat, 3.1*10**9,tf_bindaff, 1.2)
    return target_gene, alpha, beta


hg19_promoter_df = pd.read_csv("hg_promoter_regions.bed", delimiter=" ")
hg19_promoter_df = hg19_promoter_df.iloc[:, : 8]
hg19_promoter_df.columns = ['chr','start', 'end', 'name', 'score','strand','thickstart','thickend']
tf_knockout_df = pd.read_csv("raw_data/esr1_knocktf.csv")
tf_binding_df = pd.read_csv('MA0112.3.tsv', delimiter="\t", names=['chr','start','end','name', 'score', 'other', 'strand'])
target_gene = tf_knockout_df.loc[0,'TF']
tf_bindlocs_bindscores = tf_binding_df.loc[:,['chr','start', 'end','score']]

tg_arr = tf_knockout_df.iloc[0][:4].values.flatten().tolist()
tg, alpha, beta = tg_main(tg_arr, tf_bindlocs_bindscores, 100,100,100,700)

# with open('results.csv', 'w') as results:
#     results.write('tg, alpha, beta\n')
#     for i in range(len(tf_knockout_df.index)): #fix
#         tg_arr = tf_knockout_df.iloc[i][:4].values.flatten().tolist()
#         tg, alpha, beta = tg_main(tg_arr, tf_bindlocs_bindscores, 100,100,100,700)
#         results.write(f'{tg}, {alpha}, {beta}\n')


file_arr = ['','']

# return
with open('results/knock_results.csv', 'w') as file:
    file.write('tf,tg,beta')
    for tf_file in file_arr:
        tf_df = pd.read_csv(tf_file)

        for tg in tf_df[]