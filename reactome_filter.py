#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:00:29 2025

@author: harshada
"""
# reactome_filter.py
import pandas as pd
import os

def filter_reactome_interactions(reactome_file_path, gene_file_path, output_path):
    with open(gene_file_path, 'r') as f:
        gene_list = set(line.strip() for line in f if line.strip())

    df = pd.read_csv(reactome_file_path, sep='\t')
    filtered_df = df[(df["Gene1"].isin(gene_list)) & (df["Gene2"].isin(gene_list))]
    filtered_df.reset_index(drop=True, inplace=True)

    os.makedirs(output_path, exist_ok=True)
    output_file = os.path.join(output_path, "reactome_filtered_interactions.tsv")
    filtered_df.to_csv(output_file, sep='\t', index=False)
    
    return output_file
