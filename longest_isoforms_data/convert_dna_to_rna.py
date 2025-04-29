from Bio import SeqIO
import pandas as pd

# Load the refseq CSV file with longest isoforms information
refseq_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/refseq_longest_isoforms.csv'
df_refseq = pd.read_csv(refseq_file)

# Set the correct path for your FASTA file in the isoform_raw_data folder
fasta_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/isoform_raw_data/GCF_000001405.40_GRCh38.p14_genomic.fna'

# Function to extract DNA sequence from FASTA file
def extract_dna_sequence(fasta_file, chromosome, start, end):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == chromosome:
            return record.seq[start-1:end]  
    return None

# Function to convert DNA sequence to RNA
def convert_dna_to_rna(dna_sequence):
    return dna_sequence.transcribe()

# Iterate over the rows in the refseq dataframe and convert DNA to RNA
rna_sequences = []
for index, row in df_refseq.iterrows():
    chromosome = row['chromosome']
    start = row['start']
    end = row['end']
    
    # Extract DNA sequence from FASTA file
    dna_sequence = extract_dna_sequence(fasta_file, chromosome, start, end)
    
    if dna_sequence:
        # Convert the DNA sequence to RNA
        rna_sequence = convert_dna_to_rna(dna_sequence)
        rna_sequences.append(rna_sequence)
        print(f"Converted DNA to RNA for Transcript: {row['transcript_id']}")
    else:
        print(f"DNA sequence not found for Transcript: {row['transcript_id']}")

# Save RNA sequences to a new file 
rna_output_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/rna_sequences.csv'
df_refseq['rna_sequence'] = rna_sequences
df_refseq.to_csv(rna_output_file, index=False)

print(f"RNA sequences saved to {rna_output_file}")
