import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import gzip

def load_dna_fasta(fasta_path):
    sequences = {}
    with gzip.open(fasta_path, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = str(record.seq)
    return sequences

def load_refseq_csv(csv_path):
    df = pd.read_csv(csv_path)
    df['transcript_id'] = df['transcript_id'].str.strip()
    return df

def extract_transcript_sequences(dna_sequences, refseq_df):
    extracted_sequences = {}
    for _, row in refseq_df.iterrows():
        transcript_id = row['transcript_id']
        found = False
        for seq_id in dna_sequences:
            if transcript_id in seq_id:
                extracted_sequences[transcript_id] = dna_sequences[seq_id]
                found = True
                break
        if not found:
            print(f"[!] Transcript not found: {transcript_id}")
    return extracted_sequences

def dna_to_rna(dna_seq):
    return dna_seq.replace('T', 'U')

def rna_to_protein(rna_seq):
    seq_obj = Seq(rna_seq)
    return str(seq_obj.translate(to_stop=False))

def clean_protein_sequence(protein_seq):
    return protein_seq.replace("*", "")

def save_protein_sequences(protein_seqs, output_file):
    with open(output_file, 'w') as f:
        for transcript_id, protein_seq in protein_seqs.items():
            f.write(f">{transcript_id}\n{protein_seq}\n")

def full_pipeline(fasta_path, refseq_csv_path):
    print("Running sequence extraction and conversion...")

    # Step 1: Load FASTA and CSV
    dna_sequences = load_dna_fasta(fasta_path)
    print(f"DNA sequences loaded: {len(dna_sequences)}")

    refseq_df = load_refseq_csv(refseq_csv_path)
    print(f"Refseq CSV loaded: {len(refseq_df)} transcripts")

    # Step 2: Extract transcript sequences
    extracted = extract_transcript_sequences(dna_sequences, refseq_df)
    print(f"RNA sequences extracted: {len(extracted)} found.")

    # Step 3: Convert to RNA and Protein
    rna_seqs = {tid: dna_to_rna(seq) for tid, seq in extracted.items()}
    protein_seqs = {tid: rna_to_protein(rna) for tid, rna in rna_seqs.items()}
    print(f"Protein sequences converted: {len(protein_seqs)}")

    # Step 4: Clean and Save
    cleaned_proteins = {tid: clean_protein_sequence(seq) for tid, seq in protein_seqs.items()}
    print(f"Cleaned protein sequences: {len(cleaned_proteins)}")

    save_protein_sequences(cleaned_proteins, "cleaned_protein_sequences.txt")
    print("Cleaned protein sequences saved to 'cleaned_protein_sequences.txt'.")

    return cleaned_proteins

# === File Paths (update if needed) ===
if __name__ == "__main__":
    dna_fasta_file = "/Users/maryamkoddus/Documents/maryam-ko-QMUL-Msc-Project/longest_isoforms_data/isoform_raw_data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    refseq_csv_file = "/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/longest_isoforms_data/data_csv/refseq_longest_isoforms.csv"

    # Run the pipeline
    full_pipeline(dna_fasta_file, refseq_csv_file)
