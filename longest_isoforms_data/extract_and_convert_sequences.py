import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip

def extract_and_convert_to_rna(dna_fasta_file, refseq_df):
    print("Extracting and converting DNA to RNA...")
    
    # Create a dictionary for quick access to the DNA sequences by transcript ID
    dna_sequences = {}
    
    with gzip.open(dna_fasta_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            dna_sequences[record.id] = record.seq
    
    print(f"DNA sequences loaded: {len(dna_sequences)} sequences.")
    
    # Extract RNA sequences based on the longest isoforms in the refseq dataframe
    rna_sequences = []
    
    for _, row in refseq_df.iterrows():
        transcript_id = row['transcript_id']
        if transcript_id in dna_sequences:
            dna_seq = dna_sequences[transcript_id]
            rna_seq = dna_seq.transcribe()  # Convert DNA to RNA
            rna_sequences.append(str(rna_seq))
    
    print(f"RNA sequences extracted: {len(rna_sequences)} sequences found.")
    return rna_sequences

def convert_rna_to_protein(rna_sequences):
    print("Converting RNA to Protein sequences...")
    protein_sequences = []
    
    for rna in rna_sequences:
        rna_seq = Seq(rna)
        protein_seq = rna_seq.translate()  # Convert RNA to protein (amino acid sequence)
        protein_sequences.append(str(protein_seq))
    
    print(f"Protein sequences converted: {len(protein_sequences)} sequences.")
    return protein_sequences

def remove_padding(protein_sequences):
    print("Removing padding from protein sequences...")
    cleaned_protein_sequences = []
    
    for seq in protein_sequences:
        cleaned_seq = seq.replace('*', '')  # Remove any '*' characters (stop codons)
        cleaned_protein_sequences.append(cleaned_seq)
    
    print(f"Padding removed: {len(cleaned_protein_sequences)} sequences cleaned.")
    return cleaned_protein_sequences

def full_pipeline(dna_fasta_file, refseq_csv_file):
    print("Starting full pipeline...")
    refseq_df = pd.read_csv(refseq_csv_file)
    print("Refseq CSV file loaded.")
    
    # Step 1: Extract and convert to RNA
    rna_sequences = extract_and_convert_to_rna(dna_fasta_file, refseq_df)
    
    # Step 2: Convert RNA to Protein sequences
    protein_sequences = convert_rna_to_protein(rna_sequences)
    
    # Step 3: Remove padding (stop codons)
    cleaned_protein_sequences = remove_padding(protein_sequences)
    
    print("Pipeline complete.")
    return cleaned_protein_sequences

if __name__ == "__main__":
    print("Running sequence extraction and conversion...")
    
    # File paths
    dna_fasta_file = "/Users/maryamkoddus/downloads/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    refseq_csv_file = "/Users/maryamkoddus/downloads/refseq_longest_isoforms.csv"
    
    # Run the full pipeline
    cleaned_protein_sequences = full_pipeline(dna_fasta_file, refseq_csv_file)
    
    # Output the result (optionally, save to a file or print)
    print(f"Cleaned protein sequences: {len(cleaned_protein_sequences)} sequences.")
    
    # Save the cleaned protein sequences to a file
    with open("/Users/maryamkoddus/downloads/cleaned_protein_sequences.txt", "w") as output_file:
        for seq in cleaned_protein_sequences:
            output_file.write(seq + "\n")
    
    print("Cleaned protein sequences saved to 'cleaned_protein_sequences.txt'.")

