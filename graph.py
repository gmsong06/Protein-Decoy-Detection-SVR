import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from Bio.PDB import PDBParser
from residue import Residue  # Assuming your custom Residue class is defined
from pdb_reader import Protein

class Complex:
    def __init__(self, pdb_file):
        self.protein = Protein(pdb_file)
        self.res1, self.res2 = self.protein.get_interface_residues()
        self.interaction_dict = self.create_interaction_dict()

    def create_interaction_dict(self):
        interaction_dict = {}
        for i in range(len(self.res1)):
            residue_1 = Residue(self.res1[i], 1)
            residue_2 = Residue(self.res2[i], 2)
            residue_1.add_interaction(residue_2)
            residue_2.add_interaction(residue_1)
            interaction_dict[(residue_1.get_id(), residue_2.get_id())] = (residue_1, residue_2)
        return interaction_dict

    def get_interaction_dict(self):
        return self.interaction_dict

    def plot_interactions(self):
        G = nx.Graph()
        all_residues = set()

        for res1, res2 in self.interaction_dict.values():
            G.add_node(res1, label=res1.get_name(), chain=res1.chain)
            G.add_node(res2, label=res2.get_name(), chain=res2.chain)
            G.add_edge(res1, res2)
            all_residues.add(res1)
            all_residues.add(res2)

        pos = {}
        x_left = -3
        x_right = 3
        y_step = 15
        y_left = 1
        y_right = 1

        for node in G.nodes():
            if node.chain == 1:
                pos[node] = (x_left, y_left)
                y_left -= y_step
            elif node.chain == 2:
                pos[node] = (x_right, y_right)
                y_right -= y_step

        color_map = ['lightblue' if node.chain == 1 else 'red' for node in G.nodes()]

        plt.figure(figsize=(12, 9))
        nx.draw(G, pos, with_labels=False, node_size=700, node_color=color_map, edge_color="gray")
        labels = nx.get_node_attributes(G, 'label')
        nx.draw_networkx_labels(G, pos, labels, font_size=12, font_color="black")
        plt.title("2D Interactive Protein Residue Interaction Network with Chain Colors")
        plt.show()

    @staticmethod
    def calculate_rmsd(dict1, dict2):
        common_keys = set(dict1.keys()).intersection(set(dict2.keys()))
        if not common_keys:
            raise ValueError("No common interactions found to calculate RMSD.")

        rmsd_sum = 0
        for key in common_keys:
            res1_a, res1_b = dict1[key]
            res2_a, res2_b = dict2[key]

            coord1_a = np.array([atom.coord for atom in res1_a])
            coord1_b = np.array([atom.coord for atom in res1_b])
            coord2_a = np.array([atom.coord for atom in res2_a])
            coord2_b = np.array([atom.coord for atom in res2_b])

            diff_a = coord1_a - coord2_a
            diff_b = coord1_b - coord2_b

            rmsd_sum += np.sum(diff_a**2) + np.sum(diff_b**2)

        rmsd_value = np.sqrt(rmsd_sum / (len(common_keys) * 2))  # 2 for pairs
        return rmsd_value

# Example usage:
complex1 = Complex('targets/1acb_complex_H.pdb')
complex2 = Complex('targets/1acb_complex_H.pdb')

print("Interaction Dictionary for Complex 1:")
for key, value in complex1.get_interaction_dict().items():
    print(f"Interaction: {key} - Residue 1: {value[0].get_name()}, Residue 2: {value[1].get_name()}")

print("\nInteraction Dictionary for Complex 2:")
for key, value in complex2.get_interaction_dict().items():
    print(f"Interaction: {key} - Residue 1: {value[0].get_name()}, Residue 2: {value[1].get_name()}")

# Plot interactions for Complex 1
# complex1.plot_interactions()

# Calculate RMSD between the two interaction dictionaries
try:
    rmsd = Complex.calculate_rmsd(complex1.get_interaction_dict(), complex2.get_interaction_dict())
    print(f"RMSD between Complex 1 and Complex 2: {rmsd}")
except ValueError as e:
    print(e)
