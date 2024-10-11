import pandas as pd
import glob
import re

# Function to clean the gene name by removing any trailing ".numbers"
def clean_gene_name(gene_name):
    return re.sub(r'\.\d+(?:\.\d+)*$', '', gene_name)

# Get all CSV files
csv_files = glob.glob("./beta_training/GEO_expression_data/Escherichia coli/*.csv")

# Dictionary to keep track of time points and groups
time_point_groups = {}

# List to store all dataframes for merging later
dataframes = []

for csv_file in csv_files:
    # Read the CSV into a dataframe
    df = pd.read_csv(csv_file)
    
    # Clean the gene names
    df.index = df.index.map(clean_gene_name)
    
    # Process the columns (time | group_no) to manage group numbers for repeated time points
    new_columns = []
    for col in df.columns:
        time_point, group_no = col.split(" | ")
        
        # If this time point already exists, assign a new group number
        if time_point in time_point_groups:
            time_point_groups[time_point] += 1
        else:
            time_point_groups[time_point] = 1
        
        # Create new column name with updated group number
        new_columns.append(f"{time_point} | {time_point_groups[time_point]}")
    
    # Update the dataframe with the new column names
    df.columns = new_columns
    
    # Append dataframe to the list
    dataframes.append(df)

# Merge all dataframes on the gene names (rows)
combined_df = pd.concat(dataframes, axis=1)

# Replace any missing values with NaN (handled by pandas automatically)
combined_df = combined_df.where(pd.notnull(combined_df), None)

# Write the combined dataframe to a CSV
combined_df.to_csv("combined_gene_expression.csv")

print("Merged CSV created successfully!")
