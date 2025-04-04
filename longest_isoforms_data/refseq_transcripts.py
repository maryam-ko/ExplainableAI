import pandas as pd

# Load your dataset
df = pd.read_csv('longest_isoforms.csv')

# Extract the RefSeq transcript IDs (assuming they are in the second column)
refseq_ids = df['transcript_id']

# Save the IDs to a text file (one per line)
with open('refseq_transcripts.txt', 'w') as f:
    for transcript_id in refseq_ids:
        f.write(f"{transcript_id}\n")
        
print("RefSeq transcript IDs have been saved to refseq_transcripts.txt")


