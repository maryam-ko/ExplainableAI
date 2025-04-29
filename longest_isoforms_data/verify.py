import pandas as pd

# Load the GTF file
gtf_file = '/Users/maryamkoddus/Downloads/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Filter to only include transcript entries (Column 2 contains 'transcript')
df_transcripts = df_gtf[df_gtf[2] == 'transcript'].copy()

# Calculate transcript length based on start and end positions (columns 4 and 3)
df_transcripts['length'] = df_transcripts[4] - df_transcripts[3] + 1

# Extract the gene_id and transcript_id from the GTF attributes (column 8)
df_transcripts['gene_id'] = df_transcripts[8].str.extract('gene_id "([^"]+)"')
df_transcripts['transcript_id'] = df_transcripts[8].str.extract('transcript_id "([^"]+)"')

# Find the longest isoform per gene by selecting the maximum length per gene
longest_isoforms_in_gtf = df_transcripts.loc[df_transcripts.groupby('gene_id')['length'].idxmax()]

# Keep only relevant columns (gene_id, transcript_id, length)
longest_isoforms_in_gtf = longest_isoforms_in_gtf[['gene_id', 'transcript_id', 'length']]

# Load the longest isoforms from the CSV file
df_refseq = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/refseq_longest_isoforms.csv')
# Checking a few rows of the GTF file and the refseq CSV
print(df_gtf[['gene_id', 'transcript_id']].head())
print(df_refseq[['gene_id', 'transcript_id']].head())

# Merge the two dataframes on 'gene_id' and 'transcript_id'
merged_df = pd.merge(longest_isoforms_in_gtf, df_refseq, on=['gene_id', 'transcript_id'], how='outer', indicator=True)

# Show rows with discrepancies (missing in either GTF or CSV)
discrepancies = merged_df[merged_df['_merge'] != 'both']
print(discrepancies)

# Merge the dataframes to compare lengths
merged_length_check = pd.merge(longest_isoforms_in_gtf, df_refseq, on=['gene_id', 'transcript_id'], how='inner')

# Check if lengths match
length_mismatch = merged_length_check[merged_length_check['length_x'] != merged_length_check['length_y']]
print(length_mismatch)
