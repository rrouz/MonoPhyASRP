import pandas as pd

file_path = 'autophy_msa.fasta.state'

sequences = {}

with open(file_path, 'r') as file:
    current_node = None
    current_sequence = ''
    first_line = True  

    for line in file:
        if not line.startswith('#'):
            if first_line:
                first_line = False
                continue

            parts = line.strip().split('\t')
            node, state = parts[0], parts[2]

            if node != current_node:
                if current_node is not None:
                    sequences[current_node] = current_sequence  
                current_node = node
                current_sequence = state
            else:
                current_sequence += state

    if current_node is not None:
        sequences[current_node] = current_sequence

with open('ancestral_node_sequences.fasta', 'w') as fasta_file:
    for node, sequence in sequences.items():
        fasta_file.write(f'>{node}\n{sequence}\n')
