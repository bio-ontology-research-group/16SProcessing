import sys
import re

otutable = sys.argv[1]
sintax_reads = sys.argv[2]

genus_dict = {}
with open(sintax_reads, 'r') as file:
    for line in file:
        parts = line.strip().split(';')
        OTU = parts[0]
        for part in parts:
            #print(part)
            section = part.split(',')
            for sect in section:
                if sect.startswith('g:'):
                    genus = re.match(r'g:([^(\s]+)', sect)
                    genus = re.sub(r'g:', '', genus.group(0))
                    genus_dict[OTU] = genus

# Read the contents of the second file (otutab_relative.txt), replace OTU_X with the genus names,
# and write the result into a new file (otutab_relative_taxa.txt)
with open(otutable, 'r') as input_file, open('otutab_relative_withtaxa.txt', 'w') as output_file:
    header = input_file.readline()  # Read and write the header line as is
    output_file.write(header)
    
    for line in input_file:
        parts = line.strip().split('\t')
        otu_id = parts[0]
        if otu_id in genus_dict:
            genus_name = genus_dict[otu_id]
            parts[0] = genus_name  # Replace OTU_X with the genus name
            new_line = '\t'.join(parts) + '\n'
            output_file.write(new_line)
        else:
            output_file.write(line)  # If OTU_X is not found, write the line as is
