import pandas as pd
import glob
import re

# Function to clean the gene name by removing any trailing ".numbers"
def clean_gene_name(gene_name):
    return re.sub('^[a-zA-Z]+(?:\\.\\d+)$', '', gene_name)

# Get all CSV files
csv_files = glob.glob("./beta_training/GEO_expression_data/Escherichia coli/*.csv")

# Dictionary to keep track of time points and groups
time_point_groups = {}

# List to store all dataframes for merging later
dataframes = []

for csv_file in csv_files:
    # Read the CSV into a dataframe
    df = pd.read_csv(csv_file, index_col="Unnamed: 0")
    
    # Clean the gene names
    df.index = df.index.map(clean_gene_name)

    # Process the columns (time | group_no) to manage group numbers for repeated time points
    new_columns = []
    for col in df.columns:
        try:
            time_point, group_no = col.split(" | ")
        except:
            time_point, group_no = col.strip(), 0
        
        # If this time point already exists, assign a new group number
        if time_point in time_point_groups:
            time_point_groups[time_point] += 1
        else:
            time_point_groups[time_point] = 1
        
        # Create new column name with updated group number
        new_columns.append(f"{time_point} | {time_point_groups[time_point]}")
    
    # Update the dataframe with the new column names
    df.columns = new_columns

    if df.index.duplicated().any():
        duplicates = df[df.index.duplicated(keep=False)]  # Get all duplicate rows
        unique_idx = df.index.unique()  # Get all unique gene names
        
        # Create a dictionary to collect all data for each gene
        data_dict = {gene: {} for gene in unique_idx}
        
        for col in df.columns:
            time = re.search("(\\d+) | \\d+", col).group(0).strip()
            time_point_groups[time] += 1

            for gene in unique_idx:
                gene_data = df.loc[gene, col]
                
                if isinstance(gene_data, pd.Series):  # Multiple rows for this gene (duplicates)
                    data_dict[gene][f"{time} | {time_point_groups[time]}"] = gene_data.iloc[0]
                else:  # Single row for this gene
                    data_dict[gene][col] = gene_data
        
        # Convert the dictionary into a DataFrame
        resolved_df = pd.DataFrame.from_dict(data_dict, orient='index')
        
        # Append the resolved dataframe to the list of dataframes
        dataframes.append(resolved_df)
    else:
        # No duplicates, append df directly
        dataframes.append(df)

# Concatenate all dataframes into one final dataframe
combined_df = pd.concat(dataframes, axis=1)

# Write the combined dataframe to a CSV
combined_df.to_csv("combined_gene_expression.csv")

print("Merged CSV created successfully!")