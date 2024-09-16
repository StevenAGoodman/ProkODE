import maths 
import pandas as pd

# input: prob of rnap binding at instant, N_tf, Kd_tf, N_p, Kd_p, Nns
# output: beta

def main(data_file):
    input_data = pd.read_csv('results/result.csv', names=['TF_ID', 'TG_ID', 'R', 'N_tf', 'Kd_tf', 'N_p', "Kd_p", 'N_ns'])

    result_text = ''
    for index, row in input_data.iterrows:
        tf = row['TF_ID']
        tg = row['TG_ID']
        beta = maths.rev_eq1(row['R'], row['N_tf'], row['Kd_tf'], row['N_p'], row['Kd_p'], row['N_ns'])
        results += f'{tf},{tg},{beta}\n'

    
    results = open("beta_results.csv", 'w')
    results.write('tf,tg,beta')
    results.write(result_text)
    results.close()

# import GEOparse

# gse = GEOparse.get_GEO(geo="GSE1563", destdir="./")

# print()
# print("GSM example:")
# for gsm_name, gsm in gse.gsms.items():
#     print("Name: ", gsm_name)
#     print("Metadata:",)
#     for key, value in gsm.metadata.items():
#         print(" - %s : %s" % (key, ", ".join(value)))
#     print ("Table data:",)
#     print (gsm.table.head())
#     break
