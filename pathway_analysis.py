#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 13:03:51 2025

@author: harshada
"""

###############################################################################
import os
import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import textwrap

sns.set_theme(style="whitegrid")

# Load expression data
expression_file = "/home/harshada/Abhinav/Abhinav_network/results/SFA/SIM1/val/SIM1_fold_change_comparison.csv"
df_exp = pd.read_csv(expression_file)

# Extract KO and control expression
df_ko = df_exp[['genes', 'Steady State (Perturbed (GATA3)']].copy()
df_ko.columns = ['Gene', 'Expression_KO']
df_ko = df_ko[df_ko['Expression_KO'] != 0]

df_ctrl = df_exp[['genes', 'Steady State (No perturb)']].copy()
df_ctrl.columns = ['Gene', 'Expression_Control']
df_ctrl = df_ctrl[df_ctrl['Expression_Control'] != 0]

# Merge on gene names
df_merged = pd.merge(df_ko, df_ctrl, on='Gene', how='inner')
gene_list = df_merged['Gene'].tolist()

# Run enrichment
os.chdir("/home/harshada/Abhinav/Abhinav_network/results/SFA/pathway_analysis/")
enr = gp.enrichr(
    gene_list=gene_list,
    gene_sets=[
        'GO_Cellular_Component_2023',
        'CellMarker_2024',
        'GO_Molecular_Function_2023',
        'KEGG_2021_Human',
        'MSigDB_Hallmark_2020',
        'Reactome_Pathways_2024'
    ],
    organism='Human',
    outdir='enrichment_results',
    cutoff=0.05
)
df_enrich = enr.results

# Function to compute average expression
def compute_avg_expression(genes_str, col):
    genes = genes_str.split(';')
    vals = df_merged[df_merged['Gene'].isin(genes)][col]
    return vals.mean() if not vals.empty else 0

# Add columns for expression and log10 p-values
df_enrich['Avg_Expr_KO'] = df_enrich['Genes'].apply(lambda g: compute_avg_expression(g, 'Expression_KO'))
df_enrich['Avg_Expr_Control'] = df_enrich['Genes'].apply(lambda g: compute_avg_expression(g, 'Expression_Control'))
df_enrich['log10pval'] = -np.log10(df_enrich['Adjusted P-value'] + 1e-10)

# Filter by relevant biological processes
keywords = 'Signaling|Pathway|Proliferation|Cell Cycle|Migration|Invasion|EMT|TGF|PI3K|MAPK|EGFR|ERBB|mTOR|NF-κB|Wnt|Notch|BRCA|Estrogen|HER2|VEGF'
df_filtered = df_enrich[df_enrich['Term'].str.contains(keywords, case=False, regex=True)]

# Select top 25 significant pathways
top_pathways = df_filtered.sort_values(by='Adjusted P-value').head(25).copy()

# Wrap long names for readability
top_pathways['Pathway'] = top_pathways['Term'].apply(lambda x: textwrap.fill(x, width=45))

# Sort by control expression for consistent dumbbell plot
top_pathways.sort_values(by='Avg_Expr_Control', ascending=True, inplace=True)

# Plotting
fig, ax = plt.subplots(figsize=(11, 10))

# Plot dumbbell lines
for _, row in top_pathways.iterrows():
    ax.plot(
        [row['Avg_Expr_Control'], row['Avg_Expr_KO']],
        [row['Pathway'], row['Pathway']],
        color='gray', linewidth=2, alpha=0.7, zorder=1
    )

# Plot scatter points
ax.scatter(
    top_pathways['Avg_Expr_Control'], top_pathways['Pathway'],
    color='#1f77b4', edgecolor='black', s=120, label='Control', zorder=2
)
ax.scatter(
    top_pathways['Avg_Expr_KO'], top_pathways['Pathway'],
    color='#d62728', edgecolor='black', s=120, label='GATA3 KO', zorder=2
)

# Labeling and title
ax.set_xlabel("Average Expression", fontsize=14, labelpad=12)
ax.set_title("GATA3 Knockout vs Control: Expression Changes in Key Pathways", fontsize=16, pad=20)
ax.set_ylabel("")

# Styling
ax.tick_params(axis='both', labelsize=12)
sns.despine(left=True, bottom=True)
ax.grid(axis='x', linestyle='--', alpha=0.4)

# Legend
ax.legend(frameon=False, fontsize=12, loc='lower right')

# Layout & save
plt.tight_layout()
plt.savefig("expression_comparison_dumbbell_enhanced.png", dpi=700, bbox_inches='tight')
plt.show()

###############################################################################





