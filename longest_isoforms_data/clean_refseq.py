import pandas as pd

# Function to clean gene IDs (strip spaces, remove quotes)
def clean_gene_ids(df, column_name):
    df[column_name] = df[column_name].str.strip().str.replace('"', '', regex=False)
    return df

# Step 1: Clean the CSV file
csv_file = '/Users/maryamkoddus/Downloads/refseq_longest_isoforms.csv'
df_csv = pd.read_csv(csv_file)

# Clean the 'gene_id' column in the CSV file
df_csv = clean_gene_ids(df_csv, 'gene_id')

# Save the cleaned CSV
df_csv.to_csv('/Users/maryamkoddus/Downloads/cleaned_refseq_longest_isoforms.csv', index=False)
print("Gene IDs cleaned in CSV file.")

# Step 2: Clean the GTF file
gtf_file = '/Users/maryamkoddus/Downloads/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Filter only transcript entries (those with 'transcript' in column 2)
df_transcripts = df_gtf[df_gtf[2] == 'transcript'].copy()

# Extract the gene_id from the GTF attributes (column 8)
df_transcripts['gene_id'] = df_transcripts[8].str.extract('gene_id "([^"]+)"')

# Clean the 'gene_id' column in the GTF file
df_transcripts = clean_gene_ids(df_transcripts, 'gene_id')

# Save the cleaned GTF transcript data (optional)
df_transcripts.to_csv('/Users/maryamkoddus/Downloads/cleaned_GTF_transcripts.csv', index=False)
print("Gene IDs cleaned in GTF file.")

# Step 3: Check for mismatched gene IDs between the CSV and GTF files
# Load the cleaned CSV and GTF files
df_cleaned_csv = pd.read_csv('/Users/maryamkoddus/Downloads/cleaned_refseq_longest_isoforms.csv')
df_cleaned_gtf = pd.read_csv('/Users/maryamkoddus/Downloads/cleaned_GTF_transcripts.csv')

# Merge on 'gene_id' to check if the IDs match
merged_df = pd.merge(df_cleaned_csv, df_cleaned_gtf[['gene_id']], on='gene_id', how='outer', indicator=True)

# Find mismatched IDs (those present in one file but not the other)
mismatched_ids = merged_df[merged_df['_merge'] != 'both']

# Print mismatched gene IDs (if any)
if not mismatched_ids.empty:
    print(f"Mismatched gene IDs found:\n{mismatched_ids[['gene_id', '_merge']]}")
else:
    print("No mismatched gene IDs found.")

