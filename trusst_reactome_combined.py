import pandas as pd

# File paths
path = "/home/harshada/Abhinav/Abhinav_network/results/"
file1 = path + "finalized_interactions_with_regulation.tsv"
file2 = path + "trrust_final_interactions.csv"

# Load data
df1 = pd.read_csv(file1, sep='\t')  # source, target, interaction, Annotation, Regulation
df2 = pd.read_csv(file2)  # Gene1, Gene2, Annotation, Direction, Score

# Standardize columns for merging
df1 = df1.rename(columns={'source': 'Gene1', 'target': 'Gene2'})
df2 = df2.rename(columns={
    'Annotation': 'TRRUST_Annotation',
    'Direction': 'TRRUST_Direction',
    'Score': 'TRRUST_Score'
})

# Perform outer merge to include all interactions
merged_df = pd.merge(df1, df2, on=['Gene1', 'Gene2'], how='outer')

# Fill missing interaction values with TRRUST_Direction values
merged_df['interaction'] = merged_df['interaction'].fillna(merged_df['TRRUST_Direction'])

# Save the merged file
merged_file = path + "merged_interactions.tsv"
merged_df.to_csv(merged_file, sep='\t', index=False)

# # Also extract duplicates if needed
duplicates = merged_df[merged_df.duplicated(subset=['Gene1', 'Gene2'], keep=False)]
# duplicates.to_csv(path + "duplicate_interactions.tsv", sep='\t', index=False)

print(f"Merged dataframe saved to {merged_file}")
print(f"Number of duplicates: {duplicates.shape[0]}")
