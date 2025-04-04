import pandas as pd

# Load your cleaned CSV
csv_file_path = '/Users/maryamkoddus/Downloads/cleaned_refseq_longest_isoforms.csv'
df_csv = pd.read_csv(csv_file_path)

# Clean up the gene_id to match the GTF format (remove extra spaces, quotes if needed)
df_csv['gene_id'] = df_csv['gene_id'].apply(lambda x: f'gene_id "{x}"')

# Save the updated CSV
df_csv.to_csv('/Users/maryamkoddus/Downloads/refseq_longest_isoforms.csv', index=False)

print("CSV file with cleaned gene_ids has been saved!")

