import pandas as pd

# Load GTF file
gtf_file = '/Users/maryamkoddus/Downloads/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Filter for transcript entries (column 2 contains 'transcript')
df_transcripts = df_gtf[df_gtf[2] == 'transcript'].copy()

# Calculate transcript lengths (end - start + 1)
df_transcripts.loc[:, 'length'] = df_transcripts[4] - df_transcripts[3] + 1

# Extract gene_id and transcript_id from the attributes column (column 8)
df_transcripts.loc[:, 'gene_id'] = df_transcripts[8].str.extract(r'gene_id "([^"]+)"')
df_transcripts.loc[:, 'transcript_id'] = df_transcripts[8].str.extract(r'transcript_id "([^"]+)"')

# Drop rows missing gene_id or transcript_id
df_transcripts = df_transcripts.dropna(subset=['gene_id', 'transcript_id'])

# Get the longest isoform per gene (based on length)
longest_isoforms = df_transcripts.loc[df_transcripts.groupby('gene_id')['length'].idxmax()]

# Keep only relevant columns: gene_id, transcript_id, length, start, and end
longest_isoforms = longest_isoforms[['gene_id', 'transcript_id', 'length', 3, 4]]

# Rename columns for clarity
longest_isoforms.columns = ['gene_id', 'transcript_id', 'length', 'start', 'end']

# Save the result to CSV
output_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/refseq_longest_isoforms.csv'
longest_isoforms.to_csv(output_file, index=False)

print(f"Longest isoforms extracted and saved to {output_file}")















