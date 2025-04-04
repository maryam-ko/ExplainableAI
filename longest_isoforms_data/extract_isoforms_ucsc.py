import pandas as pd

# Load the UCSC table from Downloads folder (TSV format)
ucsc_file = '/Users/maryamkoddus/downloads/genome_ucsc.tsv'  # Update with actual TSV filename

# Read the TSV file, no need for skiprows since TSV has the correct header
ucsc_df = pd.read_csv(ucsc_file, sep='\t')

# Clean up column names by stripping any leading/trailing whitespace and removing special characters
ucsc_df.columns = ucsc_df.columns.str.strip().str.replace('"', '').str.replace('#', '')

# Print the columns to check if the names are fixed
print(ucsc_df.columns)

# Check the first few rows to ensure data is loaded correctly
print(ucsc_df.head())

# Ensure 'txStart' and 'txEnd' exist, otherwise adjust column names
if 'txStart' in ucsc_df.columns and 'txEnd' in ucsc_df.columns:
    # Calculate transcript length (txEnd - txStart)
    ucsc_df['length'] = ucsc_df['txEnd'] - ucsc_df['txStart']
    
    # Keep only the longest transcript per gene
    longest_ucsc = ucsc_df.loc[ucsc_df.groupby('name2')['length'].idxmax()]
    
    # Rename columns to match RefSeq format
    longest_ucsc_formatted = longest_ucsc[['name2', 'name', 'length']]
    longest_ucsc_formatted.columns = ['gene_id', 'transcript_id', 'length']
    
    # Save to CSV in Downloads folder
    output_file = '/Users/maryamkoddus/downloads/ucsc_longest_isoforms.csv'
    longest_ucsc_formatted.to_csv(output_file, index=False)
    
    print(f"Saved UCSC longest isoforms to: {output_file}")
else:
    print("Error: 'txStart' and 'txEnd' columns not found.")

