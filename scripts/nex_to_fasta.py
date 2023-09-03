import sys

input_nex_file = sys.argv[1]
output_fasta_file = sys.argv[2]

alignment_started = False
alignment_data = []

with open(input_nex_file, "r") as nex_file:
    for line in nex_file:
        line = line.strip()
        if line.startswith("MATRIX"):
            alignment_started = True
        elif alignment_started and line.startswith(";"):
            break
        elif alignment_started and line:
            parts = line.split()
            taxon = parts[0]
            sequence = "".join(parts[1:])
            alignment_data.append((taxon, sequence))

with open(output_fasta_file, "w") as fasta_file:
    for taxon, sequence in alignment_data:
        fasta_file.write(f">{taxon}\n{sequence}\n")
