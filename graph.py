import matplotlib.pyplot as plt
import networkx as nx
from pdb_reader import Protein
from residue import Residue

# Assuming your custom Protein and Residue classes are defined
# Create a Protein object
protein = Protein('targets/1mbv_complex_H.pdb')

# Get interface residues
res1, res2 = protein.get_interface_residues()

res_ids = {}
curr_ids = set()
all_residues = []

# Dictionary to store interactions with residue names
interaction_dict = {}

for i in range(len(res1)):
    residue_1 = Residue(res1[i], 1)
    residue_2 = Residue(res2[i], 2)

    residue_1.add_interaction(residue_2)
    residue_2.add_interaction(residue_1)

    if residue_1.get_id() not in curr_ids:
        all_residues.append(residue_1)
        curr_ids.add(residue_1.get_id())

    if residue_2.get_id() not in curr_ids:
        all_residues.append(residue_2)
        curr_ids.add(residue_2.get_id())

    # Add interactions to the dictionary
    interaction_dict[(residue_1.get_name(), residue_2.get_name())] = (
        residue_1, residue_2
    )

interactions = []

for residue in all_residues:
    for other in residue.get_interactions():
        if (residue, other) not in interactions and (other, residue) not in interactions:
            interactions.append((residue, other))

# Create the graph
interactions_graph = nx.Graph()

for residue in all_residues:
    interactions_graph.add_node(residue, label=residue.get_name(), chain=residue.chain)

# Use a set to avoid duplicate edges
unique_interactions = set()
for interaction in interactions:
    sorted_interaction = tuple(sorted(interaction, key=lambda x: x.get_id()))
    unique_interactions.add(sorted_interaction)

print(len(interactions))

interactions_graph.add_edges_from(interactions)

# Extract labels and chain IDs for coloring
labels = nx.get_node_attributes(interactions_graph, 'label')
chain_ids = nx.get_node_attributes(interactions_graph, 'chain')

# Manually set positions to separate nodes by chain
pos = {}
x_left = -3  # More space to the left
x_right = 3  # More space to the right
y_step = 15  # Increase y_step size for more vertical space
y_left = 1
y_right = 1

for node in interactions_graph.nodes():
    print(node.chain)
    if node.chain == 1:
        pos[node] = (x_left, y_left)
        y_left -= y_step
    elif node.chain == 2:
        pos[node] = (x_right, y_right)
        y_right -= y_step

# Create a color map based on chain ID
color_map = []
for node in interactions_graph.nodes():
    if node.chain == 1:
        color_map.append('lightblue')
    elif node.chain == 2:
        color_map.append('red')
    else:
        color_map.append('gray')  # Default color for any other chains

# Create a 2D plot
plt.figure(figsize=(12, 9))

# Plot nodes with colors based on chain ID
nx.draw(interactions_graph, pos, with_labels=False, node_size=700, node_color=color_map, edge_color="gray")

# Annotate nodes with labels
nx.draw_networkx_labels(interactions_graph, pos, labels, font_size=12, font_color="black")

# Set title
plt.title("2D Interactive Protein Residue Interaction Network with Chain Colors")

# Show plot
plt.show()

count = 1
# Print out the interaction dictionary
for key, value in interaction_dict.items():
    print(count)
    print(f"Interaction: {key} - Residue 1: {value[0].get_name()}, Residue 2: {value[1].get_name()}")
    count += 1

print(count)
