import pandas as pd

# Load both datasets (RefSeq and UCSC) into pandas DataFrames
refseq_df = pd.read_csv('/Users/maryamkoddus/downloads/refseq_longest_isoforms.csv')
ucsc_df = pd.read_csv('/Users/maryamkoddus/downloads/ucsc_longest_isoforms.csv')

# Ensure that the gene_id and transcript_id columns exist in both
print(refseq_df.head())  # Check the first few rows of the RefSeq dataset
print(ucsc_df.head())  # Check the first few rows of the UCSC dataset

# Compare the gene_ids from both datasets
common_genes = set(refseq_df['gene_id']).intersection(set(ucsc_df['gene_id']))

# Create a DataFrame to store the comparison results
comparison_results = []

for gene in common_genes:
    # Get the rows from both datasets that correspond to the gene
    refseq_gene = refseq_df[refseq_df['gene_id'] == gene]
    ucsc_gene = ucsc_df[ucsc_df['gene_id'] == gene]
    
    # Find the longest isoform in the RefSeq dataset
    longest_refseq = refseq_gene.loc[refseq_gene['length'].idxmax()]
    
    # Find the longest isoform in the UCSC dataset
    longest_ucsc = ucsc_gene.loc[ucsc_gene['length'].idxmax()]
    
    # Store the comparison result for this gene
    comparison_results.append({
        'gene_id': gene,
        'longest_refseq_transcript_id': longest_refseq['transcript_id'],
        'longest_refseq_length': longest_refseq['length'],
        'longest_ucsc_transcript_id': longest_ucsc['transcript_id'],
        'longest_ucsc_length': longest_ucsc['length']
    })

# Convert the comparison results into a DataFrame
comparison_df = pd.DataFrame(comparison_results)

# Save the comparison to a new CSV
comparison_df.to_csv('/Users/maryamkoddus/downloads/comparison_longest_isoforms.csv', index=False)

print("Comparison complete! Saved to 'comparison_longest_isoforms.csv'")

