# PhyASR
PhyASR: This pipeline aims to perform ancestral sequence reconstruction using a monophylogenetic tree generated from protein sequences.

## Pipeline Overview
The pipeline consists of several steps:

### Prepare Sequences (`replace_underscores_in_headers`)

- **Input:** protein_sequences.fasta
- **Output:** modified_protein_sequences.fasta

This step prepares protein sequence data by modifying headers for downstream analysis.

### Generate Multiple Sequence Alignment (`align_protein_sequences`)

- **Input:** modified_protein_sequences.fasta
- **Output:** aligned_protein_sequences.fasta

This step performs multiple sequence alignment using MAFFT.

### Phylogenetic Tree Inference (`iqtree`)

- **Input:** aligned_protein_sequences.fasta
- **Output:** IQ-TREE_output.treefile

This step infers a phylogenetic tree using IQ-TREE with LG+G+F model and bootstrapping.

### Convert Tree to Newick Format (`convert_to_nwk`)

- **Input:** IQ-TREE_output.treefile
- **Output:** IQ-TREE_output.nwk

This step converts the tree to Newick format.

### Activate Autophy (`activate_autophy`)

- **Input:** IQ-TREE_output.nwk

This step activates Autophy for ancestral sequence reconstruction.

### Autophy Tree Renamer (`rename_tree_file`)

- **Output:** output/autophy_tree.nwk

This step renames the Autophy output tree file to a standardized name.

### Re-label MSA Taxa (`autophy_msa`)

- **Input:** output/autophy_tree.nwk, aligned_protein_sequences.fasta
- **Output:** autophy_msa.fasta

This step re-labels taxa in the multiple sequence alignment (MSA) according to the Autophy tree.

### Ancestral Sequence Reconstruction (`iqtree_asr`)

- **Input:** output/autophy_tree.nwk, autophy_msa.fasta
- **Output:** autophy_msa.fasta.state

This step performs ancestral sequence reconstruction using IQ-TREE.

## Software Dependencies
Before running the pipeline, ensure you have the following software and dependencies installed:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/) - Multiple sequence alignment tool.
- [IQ-TREE](http://www.iqtree.org/) - Phylogenetic tree inference software.
- [Conda](https://conda.io/projects/conda/en/latest/index.html) - Package and environment management system.
- [Autophy](https://github.com/aortizsax/autophy) - Ancestral sequence reconstruction tool. Install using Conda:

  ```bash
  conda install -c conda-forge autophy
