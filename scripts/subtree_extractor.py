import re
from Bio import Entrez, SeqIO

newick_tree_file_path = "output/autophy_tree.nwk"

target_clade_number = "4"

with open(newick_tree_file_path, "r") as file:
    newick_tree_string = file.read()

taxon_pattern = r"tr_([^|]+)\|{}:".format(target_clade_number)
#taxon_pattern = r"([^,;()]*\|{}:[^,;()]+)".format(target_clade_number)

def extract_and_format_subtree(tree_string, target_number):
    matches = re.findall(taxon_pattern, tree_string)
    
    if matches:
        formatted_matches = [match.split('_')[0] for match in matches]
        return formatted_matches
    else:
        return None

selected_subtree = extract_and_format_subtree(newick_tree_string, target_clade_number)

if selected_subtree:
    print(f"Accession numbers extracted from the tree: {', '.join(selected_subtree)}")

Entrez.email = ""
db = "protein"  #"nucleotide" if you want DNA sequences

def fetch_sequences(accession_numbers):
    sequences = []

    for accession in accession_numbers:
        try:
            handle = Entrez.efetch(db=db, id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            sequences.append(record)
        except Exception as e:
            print(f"Failed to retrieve sequence for {accession}: {e}")

    return sequences

sequences = fetch_sequences(selected_subtree)

fasta_filename = f"subtree_sequences.fasta"

with open(fasta_filename, "w") as fasta_file:
    SeqIO.write(sequences, fasta_file, "fasta")

print(f"Sequences saved in '{fasta_filename}'")
