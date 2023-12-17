import pandas as pd
import sys

# Read OTU table
otu_table_file = sys.argv[1]
otu_table = pd.read_csv(otu_table_file, sep="\t")

# Rename the first column to "OTU" (assuming the first column is the OTU ID column)
otu_table.rename(columns={otu_table.columns[0]: "OTU"}, inplace=True)

# Set the new "OTU" column as index
otu_table.set_index("OTU", inplace=True)

# Calculate the sum of abundances for each sample
sample_sums = otu_table.sum(axis=0)

# Divide each cell in the table by the corresponding sample sum to get relative abundances
relative_abundance_table = otu_table.div(sample_sums, axis=1)

# Reset index to convert "OTU" back to a regular column
relative_abundance_table_reset = relative_abundance_table.reset_index()

# Save the modified relative abundance table to a new file
relative_abundance_table_reset.to_csv("otutab_relative.txt", sep="\t", index=False)
