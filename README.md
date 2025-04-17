# Network Pipeline

This repository contains a pipeline for processing gene interaction data from **Reactome** and **TRRUST** databases. It integrates gene expression data, interaction filtering, and network analysis to facilitate further simulations and analyses. Below is a description of the functionality of each script and the contents of the associated files.

## Table of Contents
1. [Reactome Codes](#reactome-codes)
2. [TRUST Codes](#trusst-codes)
3. [SFA Codes](#sfa-codes)
4. [Folder Structure](#folder-structure)

---

## Reactome Codes

1. **reactome_filter.py**  
   - Filters gene interactions for each gene with all other genes in the input list.  
   - Inputs: `filtered_log2fc_0.5_GSE72141.txt`, `reactome_interactions.txt`  
   - Output: Intermediate filtered interaction files.

2. **edge_interactions.py**  
   - Processes the edge interactions by applying activation (`->`), inhibition (`-|`), and removing predicted entries.  
   - Input: Reactome interactions from `reactome_filter.py`.  
   - Output: Processed edge interactions.

3. **reactome_network_pipeline.py**  
   - Main driver code that integrates the previous two scripts.  
   - Outputs:  
     - `finalized_interactions_with_regulation.tsv`  
     - `reactome_filtered_interactions_edges_ready.tsv`  
     - `reactome_filtered_interactions.tsv`

---

## TRRUST Codes

1. **trusst_interactions16apr.py**  
   - Filters and removes duplicates from the TRRUST interactions using `trrust_rawdata.tsv` and `filtered_log2fc_0.5_GSE72141.txt`.  
   - Outputs:  
     - `trusst_interactions_with_duplicates.csv`  
     - `trusst_final_interactions.csv`

2. **trusst_reactome_combined.py**  
   - Combines the TRRUST and Reactome interactions and checks for duplicates (none found).  
   - Inputs: `finalized_interactions_with_regulation.tsv`, `trrust_final_interactions.csv`.

---

## SFA Codes

1. **basal_log_transform_trrust_reactome.py**  
   - Transforms basal normalized counts into log values.  
   - Input: `/results/SFA/data/basal_expressions.csv`  
   - Output: `/results/SFA/data/trusst_reactome_basal_exp.csv`

2. **SFA_SIM1.py**  
   - Runs the first simulation of the Signal Flow Analysis (SFA) pipeline.  
   - Outputs: Multiple files saved in different directories under `/results/SFA/SIM1/`.

---

## Folder Structure
Refer to the documentation file for the folder structure

## Running the pipeline:

Start by running the Reactome-related code using reactome_network_pipeline.py.
Then run the TRRUST interaction filtering using trusst_interactions16apr.py.
Combine Reactome and TRRUST interactions with trusst_reactome_combined.py.
Use basal_log_transform_trrust_reactome.py to perform log transformation.
Run the SFA simulation with SFA_SIM1.py.

## Output files: All the results of the pipeline will be saved in the results/ directory. You can export the final network to Cytoscape for visualization.

## Citation
If you use this pipeline in your work, please cite the following paper:
Signal Flow Analysis (SFA): [Topological estimation of signal flow in complex signaling networks](https://www.nature.com/articles/s41598-018-23643-5)

-- This README provides a clear description of your repository's structure, the functionality of each script, and the steps to get started. 

# Thank you!
