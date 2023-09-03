import sys
from Bio import Phylo, SeqIO
import difflib

tree_file_path = sys.argv[1]
msa_file_path = sys.argv[2]
output_msa_path = sys.argv[3]

tree = Phylo.read(tree_file_path, "newick")
taxon_names = [clade.name for clade in tree.get_terminals()]

msa_records = list(SeqIO.parse(msa_file_path, "fasta"))

for record in msa_records:
    best_match = difflib.get_close_matches(record.id, taxon_names, n=1)
    if best_match:
        record.id = best_match[0]
        record.description = ""

with open(output_msa_path, "w") as output_handle:
    SeqIO.write(msa_records, output_handle, "fasta")
