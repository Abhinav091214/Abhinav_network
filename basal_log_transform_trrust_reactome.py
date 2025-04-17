#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 09:32:11 2025

@author: harshada
"""
import pandas as pd
import numpy as np
import os

raw_basal_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/data/basal_expressions.csv"
nodes_file = "/home/harshada/Abhinav/Abhinav_network/results/cytoscape/nodes_trrust_reactome"
output_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/data/"

# Read gene list from file
with open(nodes_file, "r") as f:
    gene_lines = [line.strip() for line in f if line.strip()]
    
# Load CSV file
df = pd.read_csv(raw_basal_path)
df.rename(columns = {'Unnamed: 0':'genes'},inplace=True)

# Filter rows where 'genes' column matches gene_lines
filtered_df = df[df['genes'].isin(gene_lines)].copy()  # Ensure a separate copy

# Compute mean values
filtered_df.loc[:, 'knockout_mean'] = filtered_df[['knockout_1', 'knockout_2', 'knockout_3']].mean(axis=1)
filtered_df.loc[:, 'control_mean'] = filtered_df[['wtGATA3_1','wtGATA3_2', 'wtGATA3_3']].mean(axis=1)

filtered_df["log10_knockout_mean"] = np.log10(filtered_df["knockout_mean"])
filtered_df["log10_control_mean"] = np.log10(filtered_df["control_mean"])

filtered_df.to_csv(os.path.join(output_path,"trusst_reactome_basal_exp.csv"), index=False)
