# MonoPhyASRP
Monophyletic Phylogenetic Ancestral Sequence Reconstruction Pipeline is a comprehensive bioinformatics workflow designed to perform ancestral sequence reconstruction from protein sequences. It begins by preparing protein sequence data, and conducting a multiple sequence alignment via MAFFT. A phylogenetic tree is then inferred using IQ-TREE with an LG+G+F model followed by a followed by a format converstion. Afterwhich, monophyletic clustering of the previously inferred phylogenetic tree is achieved through the use of Autophy. Finally, IQ-TREE is reemployed for ancestral sequence reconstruction. To execute this pipeline, users need to install and have proficiency with Snakemake, MAFFT, IQ-TREE, Conda, and Autophy.

## Dependencies
Before running the pipeline, ensure you have the following software installed:

- [MAFFT](https://mafft.cbrc.jp/alignment/software/) - Multiple sequence alignment tool.
- [IQ-TREE](http://www.iqtree.org/) - Phylogenetic tree inference software.
- [Conda](https://conda.io/projects/conda/en/latest/index.html) - Package and environment management system.
- [AutoPhy](https://github.com/aortizsax/autophy) - Ancestral sequence reconstruction tool.
