#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:48:23 2025

@author: harshada
"""

###############################################################################
import pandas as pd
import os

# File paths
trrust_file = "/home/harshada/Abhinav/Abhinav_network/data/trrust_rawdata.tsv"
gene_file_path = "/home/harshada/Abhinav/Abhinav_network/data/filtered_log2fc_0.5_GSE72141.txt"
output_dir = "/home/harshada/Abhinav/Abhinav_network/results/"

# Load TRRUST data
df = pd.read_csv(trrust_file, sep="\t")

# Define the headers you want to add
headers = ["Gene1", "Gene2", "Annotation", "Reference"]  # <-- update as needed

# Step 2: Assign your custom headers
df.columns = headers

df = df[df["Annotation"] != "Unknown"]

df["Direction"] = df["Annotation"].map({
    "Activation": "->",
    "Repression": "-|"
})

df["Score"] = 1

df.drop(columns="Reference",inplace=True)


# Load gene list
with open(gene_file_path, 'r') as f:
    gene_list = set(line.strip() for line in f if line.strip())


# Step 1: Filter TRRUST based on gene list
filtered_df = df[
    (df["Gene1"].isin(gene_list)) & (df["Gene1"].isin(gene_list))]

filtered_df.to_csv(os.path.join(output_dir,"trusst_interactions_with_duplicates.csv"),index=False)

###############################################################################

duplicates = filtered_df[filtered_df.duplicated(subset=["Gene1", "Gene2"], keep=False)]

###############################################################################

''' 
    Manually searched the references of the duplicates 
    and then kept only one interaction for each duplicate entry.
    
    The duplicates can be found in the TRRUST directory in a file named
    "duplicate_entries". These entries were searched in papers and only relevant interactions
    were picked which is stored in "duplicates_required".
    
    Now below the code is to extract the a filtered_df which includes only the required 
    interactions without any duplicates.
'''

import pandas as pd
import os 

# File paths
trusst_interactions_path = "/home/harshada/Abhinav/Abhinav_network/results/trusst_interactions_with_duplicates.csv"
duplicates_required_file = "/home/harshada/Abhinav/Abhinav_network/data/duplicates_required"

# Load TRRUST filtered interactions
df1 = pd.read_csv(trusst_interactions_path)

# Remove all duplicates based on Gene1 and Gene2
duplicates = df1[df1.duplicated(subset=["Gene1", "Gene2"], keep=False)]
df_final = df1[~df1.index.isin(duplicates.index)]

# Load the duplicates_required entries
required_df = pd.read_csv(duplicates_required_file, sep="\t")

# Keep only relevant columns
required_df = required_df[["Gene1", "Gene2", "Annotation"]]

# Append required duplicates to the cleaned dataframe
df_final= pd.concat([df_final, required_df], ignore_index=True)

# Add Direction column based on Annotation
df_final["Direction"] = df_final["Annotation"].map({
    "Activation": "->",
    "Repression": "-|"
})

# Add Score column
df_final["Score"] = 1


df_final.to_csv(os.path.join(output_dir,"trrust_final_interactions.csv"), index=False)

print(f"Final rows after cleaning and appending required duplicates: {len(df_final)}")

