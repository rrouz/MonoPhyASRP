def convert_header(header):
    # Assuming the header format is something like 'Q9Z9G5.1 RecName: Full=...'
    accession = header.split('.')[0]  # Extract the part before the first dot
    # Create a modified UniProt-style header
    uniprot_header = f">tr|{accession}"
    return uniprot_header

# Read the input FASTA file
input_file = "subtree_sequences.fasta"
output_file = "modified_subtree_sequences.fasta"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    header = ""
    sequence = ""
    for line in infile:
        if line.startswith(">"):
            if header:
                uniprot_header = convert_header(header)
                outfile.write(f"{uniprot_header}\n{sequence}\n")
            header = line[1:].strip()
            sequence = ""
        else:
            sequence += line.strip()
    # Write the last sequence
    if header:
        uniprot_header = convert_header(header)
        outfile.write(f"{uniprot_header}\n{sequence}\n")
