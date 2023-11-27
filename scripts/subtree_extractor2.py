import re

newick_tree_file_path = "output/autophy_tree.nwk"
target_clade_number = "4"

output_file_path = "output.txt"  # Specify the path to the output file
new_newick_tree_file_path = "new_tree.nwk"  # Specify the path for the new Newick file

with open(newick_tree_file_path, "r") as file:
    newick_tree_string = file.read()

taxon_pattern = r"([^,;()]*\|{}:[^,;()]+)".format(target_clade_number)

def extract_and_format_subtree(tree_string, target_number):
    matches = re.findall(taxon_pattern, tree_string)

    if matches:
        formatted_matches = []

        # Find the index of the first occurrence of the target clade
        first_occurrence_index = tree_string.find(matches[0])

        if first_occurrence_index >= 0:
            # Traverse backward to find the opening parenthesis
            while first_occurrence_index > 0 and tree_string[first_occurrence_index] != '(':
                first_occurrence_index -= 1
            
            if first_occurrence_index > 0:
                # Traverse forward to count parentheses and extract the subtree
                open_parentheses_count = 0
                subtree = ""
                for char in tree_string[first_occurrence_index:]:
                    subtree += char
                    if char == '(':
                        open_parentheses_count += 1
                    elif char == ')':
                        open_parentheses_count -= 1
                        if open_parentheses_count == 0:
                            break

                formatted_subtree = subtree

                return formatted_subtree

    return None

selected_subtree = extract_and_format_subtree(newick_tree_string, target_clade_number)

if selected_subtree:
    
    with open(new_newick_tree_file_path, "w") as newick_file:
        newick_file.write(selected_subtree + ");")

    print(selected_subtree + ");")
else:
    print(f"No clades with the specified target number '{target_clade_number}' found in the tree.")
