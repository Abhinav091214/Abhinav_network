#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 23 11:38:42 2025

@author: harshada
"""
import os
import pandas as pd
from collections import defaultdict

def load_edges(file_path):
    try:
        df = pd.read_csv(file_path, delimiter=',')
        edges = defaultdict(list)
        edge_signals = {}

        for _, row in df.iterrows():
            gene1 = row["Source"].strip().upper()
            gene2 = row["Target"].strip().upper()
            edges[gene1].append(gene2)
            edge_signals[(gene1, gene2)] = row["Signal Flow"]

        return edges, edge_signals

    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return defaultdict(list), {}

def build_edges_tree(gene, edges, visited=None):
    if visited is None:
        visited = set()
    gene = gene.upper()
    visited.add(gene)
    edges_tree = {}

    if gene in edges:
        for partner in edges[gene]:
            if partner not in visited:
                edges_tree[partner] = build_edges_tree(partner, edges, visited)

    return edges_tree

def flatten_tree(tree, edge_signals, parent_path=[]):
    flat_data = []

    for key, subtree in tree.items():
        current_path = parent_path + [key]
        flat_data.append((current_path, edge_signals.get((parent_path[-1], key)) if parent_path else None))
        flat_data.extend(flatten_tree(subtree, edge_signals, current_path))

    return flat_data

def save_edges_tree_to_excel(tree, edge_signals, file_path):
    flat_data = flatten_tree(tree, edge_signals)

    if not flat_data:
        print("No edges found.")
        return

    max_depth = max(len(path) for path, _ in flat_data)
    df = pd.DataFrame([path + [signal] for path, signal in flat_data],
                      columns=[f'Level {i+1}' for i in range(max_depth)] + ['Signal Flow'])
    df.to_excel(file_path, index=False)
    print(f"Edges tree saved to {file_path}")

# === FILE INPUTS ===
file_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/signal_flow/SIM1_signal_flow_inhibition.csv"
edges, edge_signals = load_edges(file_path)

# genes_df = pd.read_csv("/data/SFA/SFA_microgravity/Completed_sims/Receptor_identification/Candidate_receptors/SFA_Data_2024-06-08_19-03-29.964426/Perturbed_Genes.csv")
# genes = genes_df.iloc[:, 0].dropna().tolist()

genes = ["GATA3"]
output_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/path_traversal/GATA3_path/"

for gene in genes:
    edges_tree = build_edges_tree(gene, edges)

    if edges_tree:
        save_edges_tree_to_excel(edges_tree, edge_signals, os.path.join(output_path,f'{gene}_Edges_tree.xlsx'))
    else:
        print(f"No edges found for {gene}")
        