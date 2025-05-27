#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 15 10:07:57 2025

@author: harshada
"""
import pandas as pd
import powerlaw
import matplotlib.pyplot as plt
import numpy as np

# Load the data
df = pd.read_csv("/home/harshada/Abhinav/Abhinav_network/results/cytoscape/centrality.csv")
degree_data = df['Degree']

# Fit a power-law model (discrete data)
fit = powerlaw.Fit(degree_data, discrete=True)
print(f"Estimated xmin: {fit.xmin}")
print(f"Estimated alpha: {fit.power_law.alpha:.4f}")

# Extract empirical CCDF
x_empirical, y_empirical = fit.ccdf()
x_empirical = np.array(x_empirical)
y_empirical = np.array(y_empirical)

# Select ~25 log-spaced indices for sampled points
indices = np.unique(np.round(np.logspace(0, np.log10(len(x_empirical) - 1), 25)).astype(int))
x_sampled = x_empirical[indices]
y_sampled = y_empirical[indices]

# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot empirical data (25 points)
ax.scatter(x_sampled, y_sampled, color='black', s=50, label='Network Data', zorder=3)

# Plot power-law fit
fit.power_law.plot_ccdf(ax=ax, color='blue', linestyle='--', linewidth=2, label='Power-law Fit')

# Log-log scaling
ax.set_xscale('log')
ax.set_yscale('log')

# Axis labels
ax.set_xlabel("Degree", fontsize=16, fontweight='bold')
ax.set_ylabel("Probability Density", fontsize=16, fontweight='bold')

# Title
# ax.set_title("Degree Distribution with Power-law Fit", fontsize=16, fontweight='bold')

# Axis tick formatting
ax.tick_params(axis='both', which='major', labelsize=14, width=1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(1.5)
ax.spines['bottom'].set_linewidth(1.5)

# Grid and legend
ax.grid(True, which='major', linestyle='--', linewidth=0.5, alpha=0.7)
ax.legend(fontsize=12, frameon=False, loc='best')

# Save the figure
plt.tight_layout()
plt.savefig("/home/harshada/Abhinav/Abhinav_network/results/cytoscape/degree_distribution_powerlaw_enhanced.png", dpi=900, bbox_inches='tight')
plt.show()


##############################################################################
import os
import pandas as pd

edge_table = "/home/harshada/Abhinav/Abhinav_network/results/cytoscape/main_paths_edge_table.csv"

df = pd.read_csv(edge_table)

# Extract source, interaction, target
df[['Source', 'Interaction', 'Target']] = df['name'].str.extract(r'^(.*?) \((.*?)\) (.*)$')

# Add Type and Weight
df['Type'] = 'Directed'
df['Weight'] = 1

# Keep only required columns
df = df[['Source', 'Target', 'Type', 'Interaction', 'Weight']]

df.to_csv("/home/harshada/Abhinav/Abhinav_network/results/cytoscape/main_paths_source_trg.csv")


# Combine unique Source and Target nodes
nodes = pd.unique(df[['Source', 'Target']].values.ravel())

# Create node table
node_table = pd.DataFrame({'Id': nodes, 'Label': nodes})

# Export to CSV
node_table.to_csv("/home/harshada/Abhinav/Abhinav_network/results/cytoscape/gephi_node_table.csv", index=False)

################################################################################