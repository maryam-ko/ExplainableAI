import pandas as pd
from Bio import SeqIO

# Load GTF file and extract chromosome information
gtf_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/isoform_rawdata/GCF_000001405.40_GRCh38.p14_genomic.gtf'
df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

# Filter for transcript entries and extract chromosome (seqname)
df_transcripts = df_gtf[df_gtf[2] == 'transcript'].copy()
df_transcripts['chromosome'] = df_transcripts[0]  # Column 0 contains the chromosome names (seqname)

# Extract relevant columns: gene_id, transcript_id, start, end, chromosome
df_transcripts['gene_id'] = df_transcripts[8].str.extract(r'gene_id "([^"]+)"')
df_transcripts['transcript_id'] = df_transcripts[8].str.extract(r'transcript_id "([^"]+)"')
df_transcripts = df_transcripts[['gene_id', 'transcript_id', 3, 4, 'chromosome']]  # 3=Start, 4=End

# Calculate length of each transcript
df_transcripts['length'] = df_transcripts[4] - df_transcripts[3] + 1

# Get the longest isoform for each gene
longest_isoforms = df_transcripts.loc[df_transcripts.groupby('gene_id')['length'].idxmax()]

# Save the longest isoforms to a CSV file
output_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/refseq_longest_isoforms.csv'
longest_isoforms = longest_isoforms[['gene_id', 'transcript_id', 'length', 3, 4, 'chromosome']]
longest_isoforms.columns = ['gene_id', 'transcript_id', 'length', 'start', 'end', 'chromosome']
longest_isoforms.to_csv(output_file, index=False)

print(f"Longest isoforms with chromosome extracted and saved to {output_file}")















