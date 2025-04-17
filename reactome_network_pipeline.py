#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 18:14:57 2025

@author: harshada
"""
# pipeline.py
import os
from reactome_filter import filter_reactome_interactions
from edge_interactions import process_edge_directions
import pandas as pd

# === Working directory ===
base_dir = "/home/harshada/Abhinav/Abhinav_network/"

# === Paths ===
reactome_file = os.path.join(base_dir, "data", "reactome_interactions.txt")
gene_file = os.path.join(base_dir, "data", "filtered_log2fc_0.5_GSE72141.txt")
results_dir = os.path.join(base_dir, "results")

# === Step 1: Filter reactome interactions ===
filtered_file = filter_reactome_interactions(reactome_file, gene_file, results_dir)
print(f"✅ Filtered interactions saved to: {filtered_file}")

# === Step 2: Process edge interactions- This is an intermediate file ===
df = pd.read_csv(filtered_file, sep="\t")
df.rename(columns={"Gene1": "source", "Gene2": "target"}, inplace=True)
edge_ready_file = os.path.join(results_dir, "reactome_filtered_interactions_edges_ready.tsv")
df.to_csv(edge_ready_file, sep="\t", index=False)

final_output = process_edge_directions(edge_ready_file, results_dir)
print(f"✅ Final interactions with regulation saved to: {final_output}")
###############################################################################
