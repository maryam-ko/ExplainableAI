import pandas as pd

# Load the GTF file
gtf_file = '/Users/maryamkoddus/Downloads/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'

# Read the GTF file, skipping the header comments
df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Filter to only include 'transcript' entries (column 2 contains 'transcript')
df_transcripts = df_gtf[df_gtf[2] == 'transcript'].copy()

# Calculate transcript length based on the start and end positions (columns 4 and 3)
df_transcripts.loc[:, 'length'] = df_transcripts[4] - df_transcripts[3] + 1

# Extract the gene_id and transcript_id from the GTF attributes (column 8)
df_transcripts.loc[:, 'gene_id'] = df_transcripts[8].str.extract('gene_id "([^"]+)"')
df_transcripts.loc[:, 'transcript_id'] = df_transcripts[8].str.extract('transcript_id "([^"]+)"')

# Remove rows where transcript_id is missing
df_transcripts = df_transcripts.dropna(subset=['transcript_id'])

# Find the longest transcript per gene by selecting the maximum length per gene
longest_isoforms = df_transcripts.loc[df_transcripts.groupby('gene_id')['length'].idxmax()]

# Keep only relevant columns (gene_id, transcript_id, length)
longest_isoforms = longest_isoforms[['gene_id', 'transcript_id', 'length']]

# Save the output to a CSV file
longest_isoforms.to_csv('refseq_longest_isoforms.csv', index=False)

print("Longest isoforms extracted and saved to refseq_longest_isoforms.csv")

















