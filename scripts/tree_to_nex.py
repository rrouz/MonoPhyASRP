import sys

input_tree_file = sys.argv[1]
output_nex_file = sys.argv[2]

with open(input_tree_file, "r") as tree_file:
    tree_data = tree_file.read().strip()

# Process tree_data to generate NEXUS format content
nex_content = f"# NEXUS\nbegin trees;\n\ttree tree_1 = {tree_data};\nend;"

with open(output_nex_file, "w") as nex_file:
    nex_file.write(nex_content)
