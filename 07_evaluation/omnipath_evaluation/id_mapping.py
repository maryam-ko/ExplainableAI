import pandas as pd
from Bio import SeqIO

# === Step 1: Parse UniProt FASTA to build UniProt ID → Gene Symbol mapping ===

def build_uniprot_mapping(fasta_path):
    mapping = {}
    with open(fasta_path, 'r') as fasta_file:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            header = record.description
            try:
                uniprot_id = header.split('|')[1]
                gene_parts = [part for part in header.split() if part.startswith('GN=')]
                if gene_parts:
                    gene_name = gene_parts[0].split('=')[1]
                    mapping[uniprot_id] = gene_name
            except IndexError:
                continue
    return mapping

# === Step 2: Load reference interaction file and map gene names ===

def annotate_reference(reference_path, fasta_path, output_path=None):
    print("Building UniProt ↔ Gene Symbol mapping...")
    mapping = build_uniprot_mapping(fasta_path)

    print("Loading reference interactions...")
    df = pd.read_csv(reference_path, sep='\t')  # Change sep if needed

    print("Mapping gene names to UniProt IDs...")
    df['source_gene_name'] = df['source'].map(mapping)
    df['target_gene_name'] = df['target'].map(mapping)

    print("Done! Sample with new columns:")
    print(df[['source', 'source_gene_name', 'target', 'target_gene_name']].head())

    if output_path:
        df.to_csv(output_path, sep='\t', index=False)
        print(f"Annotated file saved to {output_path}")

    return df

# === Run ===

# Replace with your actual file paths
reference_file = "/Users/maryamkoddus/Downloads/ExplainableAI/07_evaluation/omnipath_evaluation/omnipath_interactions.txt"
fasta_file = "/Users/maryamkoddus/Downloads/ExplainableAI/01_input_data/raw_data/UP000005640_9606.fasta"
output_file = "/Users/maryamkoddus/Downloads/ExplainableAI/07_evaluation/omnipath_evaluation/omnipath_human_interactions.csv"

annotate_reference(reference_file, fasta_file, output_file)