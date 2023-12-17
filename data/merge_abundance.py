import sys
import pandas as pd

# Load the data from the input file
input_file = sys.argv[1]
df = pd.read_csv(input_file, delimiter="\t")

# Group the data by the 'Genus' column and sum the abundance for each sample
grouped_df = df.groupby('OTU').sum()

# Reset the index to have 'Genus' as a column
grouped_df.reset_index(inplace=True)

# Create the output file name
output_file = 'otutab_relative_withtaxa_merged.tsv'

# Save the merged data to the output file
grouped_df.to_csv(output_file, sep='\t', index=False)