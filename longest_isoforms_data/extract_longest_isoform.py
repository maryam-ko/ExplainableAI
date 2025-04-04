import pandas as pd

#Load GTF file
gtf_file = '/Users/maryamkoddus/Downloads/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz'
df_gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None)

#Filter for transcripts 
df_transcripts = df_gtf[df_gtf[2] == 'transcript'].copy()

#Calculate transcript lengths 
df_transcripts.loc[:, 'length'] = df_transcripts[4] - df_transcripts[3] + 1

#Extract gene_id and transcript_id
df_transcripts.loc[:, 'gene_id'] = df_transcripts[8].str.extract(r'gene_id "([^"]+)"')
df_transcripts.loc[:, 'transcript_id'] = df_transcripts[8].str.extract(r'transcript_id "([^"]+)"')

#Drop rows missing IDs 
df_transcripts = df_transcripts.dropna(subset=['gene_id', 'transcript_id'])

#Get longest isoform per gene 
longest_isoforms = df_transcripts.loc[df_transcripts.groupby('gene_id')['length'].idxmax()]

#Keep only needed columns 
longest_isoforms = longest_isoforms[['gene_id', 'transcript_id', 'length']]

# Save to CSV 
output_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/refseq_longest_isoforms.csv'
longest_isoforms.to_csv(output_file, index=False)

print(f"Longest isoforms extracted and saved to {output_file}")

















