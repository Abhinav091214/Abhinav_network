#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 17 09:45:15 2025

@author: harshada
"""
############# Creating an algorithm object#####################################
import os
import sfa
import numpy as np
import networkx as nx
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from warnings import filterwarnings
filterwarnings(action='ignore', category=DeprecationWarning, message='`np.bool` is a deprecated alias')

class network(sfa.base.Data):
    def __init__(self):
        super().__init__()
        self._abbr = "main_network"
        self._name = "extracted_network"
        
        # Specify the file path for network file in SIF.
        file_name = "/home/harshada/Abhinav/Abhinav_network/results/cytoscape/main_reactome_trrust.sif"

        # Use read_sif function.
        A, n2i, dg = sfa.read_sif(file_name, as_nx=True)
        self._A = A
        self._n2i = n2i
        self._dg = dg
        self._inputs = {}
        self._i2n = {idx: name for name, idx in n2i.items()}
    # end of def __init__
# end of def class

data = network()
algs = sfa.AlgorithmSet()
alg = algs.create('SP')
data.n2i

############## Analyze the data with the algorithm ######################
alg.params.alpha = 0.5
alg.params.apply_weight_norm = True
alg.data = data  # Assign the data object to the algorithm.
alg.initialize()  # Initialize the algorithm object.
###############################################################################

################ Load Expression Data and Assign xs1 #####################
basal_exp_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/data/trusst_reactome_basal_exp.csv"

# Read the CSV file
df = pd.read_csv(basal_exp_path)

# Initialize perturbation array b with zeros
b = np.zeros((data.dg.number_of_nodes(),), dtype=float)

# Initialize xs1 with zeros (same size as b)
xs1_no_perturb = np.zeros_like(b)

# Assign log2FoldChange values to xs1
for _, row in df.iterrows():
    gene = row['genes']
    if gene in data.n2i:  # Check if the gene exists in the network
        xs1_no_perturb[data.n2i[gene]] = row['log10_control_mean']
        
###############################################################################
steady_state_file_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/"

# Define the main folder and subfolders
main_folder = os.path.join(steady_state_file_path, "SIM1")
subfolders = ["nodal_vals", "signal_flow", "val"]

# Create the main folder
os.makedirs(main_folder, exist_ok=True)

# Create subfolders inside "simulation folder (sim1 etc)"
for subfolder in subfolders:
    os.makedirs(os.path.join(main_folder, subfolder), exist_ok=True)
    
print(f"Folders created inside: {main_folder}")

################ COMPUTE STEADY STATE VALUES ##################################

xs1_GATA3_no_perturb = alg.compute(xs1_no_perturb) + alg.compute(b)
xs1_GATA3_no_perturb_df = pd.DataFrame({'Node': list(data.i2n.values()), 'Steady State (No perturb)': xs1_GATA3_no_perturb})

nodal_vals_path_noperturb = os.path.join(main_folder, "nodal_vals", "SIM1_xs1_noperturb.csv") #ALWAYS EDIT QUOTES
xs1_GATA3_no_perturb_df.to_csv(nodal_vals_path_noperturb, index=False)

# b[data.n2i['GATA3']] = 0.67 # Activation perturbation

## OR

b[data.n2i['GATA3']] = -1 # Inhibition perturbation

node_exp = xs1_GATA3_no_perturb + b

xs1_GATA3_perturbed = alg.compute(node_exp)
xs1_GATA3_perturbed_df = pd.DataFrame({'Node': list(data.i2n.values()), 'Steady State (Perturbed (GATA3)': xs1_GATA3_perturbed})

nodal_vals_path_perturb = os.path.join(main_folder, "nodal_vals", "SIM1_xs1_inhibtion.csv") #ALWAYS EDIT QUOTES
xs1_GATA3_perturbed_df.to_csv(nodal_vals_path_perturb, index=False)

###############################################################################

results_GATA3 = [] 

W1 = alg.W.copy()  # Get a copy of the weight matrix (same for all conditions)

# Compute signal flows for each basal condition
F1 = W1 * xs1_GATA3_perturbed

# # Collecting signal flows for each basal condition
ir, ic = data.A.nonzero()

for i in range(ir.size):
    idx_trg, idx_src = ir[i], ic[i]
    src = data.i2n[idx_src]
    trg = data.i2n[idx_trg]
    sf = F1[idx_trg, idx_src]  # Signal flow for the current basal condition
    results_GATA3.append({'Source': src, 'Target': trg, 'Signal Flow': sf})
    
ind_up_xs1 = np.where(xs1_GATA3_perturbed > 0)[0]  # Indices of upregulated nodes
ind_dn_xs1 = np.where(xs1_GATA3_perturbed < 0)[0]  # Indices of downregulated nodes

############ Save signal flows as a CSV file ##################################

signal_flow_file_path = os.path.join(main_folder, "signal_flow", "SIM1_signal_flow_inhibition.csv") # ALWAYS EDIT QUOTES 
NETsf_result_activation_GATA3 =  pd.DataFrame(results_GATA3)
NETsf_result_activation_GATA3.to_csv(signal_flow_file_path, index=False)

################################################################################



################# MODEL VALIDATIONS ##########################################

import pandas as pd

# Load the files
file1_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/nodal_vals/SIM1_xs1_inhibtion.csv"
file2_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/data/trusst_reactome_basal_exp.csv"

df1 = pd.read_csv(file1_path)  # Adjust delimiter if necessary
df2 = pd.read_csv(file2_path)

# Rename columns for consistency
df1.rename(columns={"Node": "genes", "Steady State (Perturbed (GATA3)": "Steady_State"}, inplace=True)

# Merge the dataframes on the "genes" column
merged_df = pd.merge(df1, df2, on="genes", how="inner")

# Create a column to compare the signs
merged_df["Sign_Match"] = (merged_df["Steady_State"] * merged_df["log10_knockout_mean"]) > 0

# Count total matched and unmatched signs
total_comparisons = len(merged_df)
matched_signs = merged_df["Sign_Match"].sum()
unmatched_signs = total_comparisons - matched_signs

# Identify positive and negative cases
positive_cases = merged_df["log10_knockout_mean"] > 0

negative_cases = merged_df["log10_knockout_mean"] < 0

# Count positive and negative matches
positive_matches = ((merged_df["Steady_State"] > 0) & (merged_df["log10_knockout_mean"] > 0)).sum()
negative_matches = ((merged_df["Steady_State"] < 0) & (merged_df["log10_knockout_mean"] < 0)).sum()

# Calculate percentages
matched_percent = (matched_signs / total_comparisons) * 100 if total_comparisons > 0 else 0
positive_match_percent = (positive_matches / positive_cases.sum()) * 100 if positive_cases.sum() > 0 else 0
negative_match_percent = (negative_matches / negative_cases.sum()) * 100 if negative_cases.sum() > 0 else 0

# Print results
print(f"Total comparisons: {total_comparisons}")
print(f"Matched signs: {matched_signs} ({matched_percent:.2f}%)")
print(f"Unmatched signs: {unmatched_signs}")
print(f"Positive sign matches: {positive_matches}/{positive_cases.sum()} ({positive_match_percent:.2f}%)")
print(f"Negative sign matches: {negative_matches}/{negative_cases.sum()} ({negative_match_percent:.2f}%)")

# Save the result
output_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/val/SIM1_model_val.csv"
merged_df.to_csv(output_path, index=False)

print(f"File saved at: {output_path}")
# ###############################################################################

########## FOLD CHANGE COMPARISON #############################################
import os
import pandas as pd

no_perturb_file = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/nodal_vals/SIM1_xs1_noperturb.csv"
perturb_file = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/nodal_vals/SIM1_xs1_inhibtion.csv"
output_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/nodal_vals/SIM1_fold_change.csv"

# Load the CSV files
no_perturb_df = pd.read_csv(no_perturb_file)
perturb_df= pd.read_csv(perturb_file)

# Merge on 'genes' while keeping all distinct columns
fold_change_merged_df = no_perturb_df.merge(perturb_df, on='Node', suffixes=('_no_perturb', '_perturb'))

# Calculate SFA_fold_change
fold_change_merged_df['SFA_fold_change'] = fold_change_merged_df['Steady State (Perturbed (GATA3)'] - fold_change_merged_df['Steady State (No perturb)'] 
fold_change_merged_df.to_csv(output_path,index=False)
###############################################################################

import pandas as pd 

# File paths
SFA_fold_change = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/nodal_vals/SIM1_fold_change.csv"
rna_seq_fold_change = "/home/harshada/Abhinav/Abhinav_network/results/SFA/data/new_GSE72141_DESEQ2_results.csv"

# Load datasets
sfa_fold_change_df = pd.read_csv(SFA_fold_change)
rna_seq_fold_change_df = pd.read_csv(rna_seq_fold_change)

# Rename columns for consistency
sfa_fold_change_df.rename(columns={"Node": "genes"}, inplace=True)

# Merge dataframes on "genes" and keep only relevant columns
merged_df1 = pd.merge(sfa_fold_change_df, rna_seq_fold_change_df[['genes', 'log2FoldChange']], on="genes", how="inner")

# Create a column to compare the signs
merged_df1["Sign_Match"] = (merged_df1['SFA_fold_change'] * merged_df1['log2FoldChange']) > 0

# Count total matched and unmatched signs
total_comparisons = len(merged_df1)
matched_signs = merged_df1["Sign_Match"].sum()
unmatched_signs = total_comparisons - matched_signs

# Identify positive and negative cases
positive_cases =  merged_df1['log2FoldChange'] > 0

negative_cases = merged_df1['log2FoldChange'] < 0

# Count positive and negative matches
positive_matches = ((merged_df1['SFA_fold_change'] > 0) & (merged_df1['log2FoldChange'] > 0)).sum()
negative_matches = ((merged_df1['SFA_fold_change'] < 0) & (merged_df1['log2FoldChange'] < 0)).sum()

# Calculate percentages
matched_percent = (matched_signs / total_comparisons) * 100 if total_comparisons > 0 else 0
positive_match_percent = (positive_matches / positive_cases.sum()) * 100 if positive_cases.sum() > 0 else 0
negative_match_percent = (negative_matches / negative_cases.sum()) * 100 if negative_cases.sum() > 0 else 0

# Print results
print(f"Total comparisons: {total_comparisons}")
print(f"Matched signs: {matched_signs} ({matched_percent:.2f}%)")
print(f"Unmatched signs: {unmatched_signs}")
print(f"Positive sign matches: {positive_matches}/{positive_cases.sum()} ({positive_match_percent:.2f}%)")
print(f"Negative sign matches: {negative_matches}/{negative_cases.sum()} ({negative_match_percent:.2f}%)")

# Save the result
output_path = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/val/SIM1_fold_change_comparison.csv"
merged_df1.to_csv(output_path, index=False)



