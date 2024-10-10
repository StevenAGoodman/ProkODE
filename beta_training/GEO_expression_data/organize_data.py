import sys
import os
import json
import re
import pandas as pd

data_dir = "/workspaces/prokode/beta_training/GEO_expression_data/"

files_list = []
for file_path in os.listdir(data_dir):
    try:
        files_list.extend([file_path + '/' + file for file in os.listdir(data_dir + "/" + file_path)])
    except:
        files_list.append(file_path)


files_list = [ data_file for data_file in files_list if re.search("^.*.csv$", data_file) != None]
print(files_list)

prev_end_val = 0

big_df = pd.DataFrame([])

for data_file in files_list[2:]:
    print("yay")
    # get data
    data_df = pd.read_csv(f"{data_dir}/{data_file}")
    print('yay2')
    data_df = data_df.rename(columns={list(data_df.columns)[0]:"geneid"})
    gene_col = data_df["geneid"]

    # normalize as concentations
    data_df = data_df.iloc[:,1:].multiply(1 / 6.02e23)

    # filter out groups
    for i in range(len(data_df.axes[1])):
        colname = data_df.columns[i]
        group_n = re.search("\\| (\\d+)", colname)
        if group_n == None:
            group_n = "0"
            replace_n = str(int(group_n) + prev_end_val)

            data_df = data_df.rename(columns={colname: colname + ' | ' + replace_n})
        else:
            group_n = group_n.group(1)
        
            replace_n = str(int(group_n) + prev_end_val)
            data_df = data_df.rename(columns = {colname: colname.replace(group_n, replace_n)})

    prev_end_val = int(replace_n)

    search_df = data_df.loc[[ re.search("\\.\\d+", gene) != None for gene in list(gene_col)],:]
    if len(search_df) != 0:
        def plus_one(x):
            id = int(re.search(" \\| (\\d+)$", x).group(1))
            x = re.sub("\\.\\d+", "", x)
            return x.replace(str(id), str(id + 1))
        search_df = search_df.rename(plus_one, axis="columns" )
        data_df = data_df.loc[[re.search("\\.\\d+", gene) == None for gene in list(gene_col)]]
        prev_end_val += 1
        print("yay???")
        new_gene_col = gene_col[[re.search("\\.\\d+", gene) != None for gene in list(gene_col)]]
        new_gene_col = [gene[:gene.find(".")] for gene in list(new_gene_col)]
        search_df.insert(0,"geneid",  new_gene_col, allow_duplicates=True)
        if str(type(big_df)) == "<class 'NoneType'>" or len(big_df) == 0:
            big_df = search_df
        else:
            big_df = big_df.merge(search_df, "outer", "geneid")


    data_df.insert(0,"geneid",  pd.Series(gene_col), allow_duplicates=True)

    if str(type(big_df)) == "<class 'NoneType'>" or len(big_df) == 0:
        big_df = data_df
    else:
        big_df = big_df.merge(data_df, "outer", "geneid")

big_df.drop_duplicates()

big_df.to_csv("/workspaces/prokode/beta_training/GEO_expression_data/data.csv")