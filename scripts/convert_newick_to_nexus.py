from ete3 import Tree

input_newick_file = snakemake.input[0]  #Get the input Newick file from Snakemake
output_nexus_file = snakemake.output[0]  #Get the output Nexus file from Snakemake

with open(input_newick_file, "r") as f:
    newick_tree_str = f.read()

newick_tree = Tree(newick_tree_str)
newick_tree.write(format=1, outfile=output_nexus_file)
