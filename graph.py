import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from Bio.PDB import PDBParser
from Bio.PDB.vectors import rotaxis2m, Vector

class Residue:
    # Assuming this is a placeholder for your custom Residue class
    def __init__(self, residue, chain):
        self.residue = residue
        self.chain = chain
        self.interactions = []

    def get_name(self):
        return self.residue.get_resname()

    def get_id(self):
        return self.residue.get_id()[1]

    def get_atoms(self):
        return self.residue.get_atoms()

    def add_interaction(self, other_residue):
        self.interactions.append(other_residue)

    def get_interactions(self):
        return self.interactions

class Protein:
    def __init__(self, pdb_file):
        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure('protein', pdb_file)

    def get_interface_residues(self):
        coordinates = []
        chains = []
        residues = []

        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    if residue.id[0] == ' ':
                        for atom in residue:
                            chains.append(chain_id)
                            coordinates.append(atom.coord)
                            residues.append(residue)

        if not chains:
            return [], []

        first_chain = chains[0]
        listA, listB = [], []
        residues_A, residues_B = [], []

        for i, chain in enumerate(chains):
            if chain == first_chain:
                listA.append(coordinates[i])
                residues_A.append(residues[i])
            else:
                listB.append(coordinates[i])
                residues_B.append(residues[i])

        dist_thresh = 5
        res_in_contact_A = []
        res_in_contact_B = []

        for i, posA in enumerate(listA):
            for j, posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                    res_in_contact_A.append(residues_A[i])
                    res_in_contact_B.append(residues_B[j])

        return res_in_contact_A, res_in_contact_B

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
    def kabsch_rmsd(P, Q):
        """
        Rotate matrix P unto matrix Q using Kabsch algorithm and calculate RMSD.
        """
        C = np.dot(np.transpose(P), Q)
        V, S, W = np.linalg.svd(C)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        U = np.dot(V, W)
        P = np.dot(P, U)

        return np.sqrt(np.sum((P - Q)**2) / len(P))

    @staticmethod
    def calculate_rmsd(dict1, dict2):
        common_keys = set(dict1.keys()).intersection(set(dict2.keys()))
        if not common_keys:
            raise ValueError("No common interactions found to calculate RMSD.")

        rmsd_sum = 0
        for key in common_keys:
            res1_a, res1_b = dict1[key]
            res2_a, res2_b = dict2[key]

            coords1_a = np.array([atom.coord for atom in res1_a.residue.get_atoms()])
            coords1_b = np.array([atom.coord for atom in res1_b.residue.get_atoms()])
            coords2_a = np.array([atom.coord for atom in res2_a.residue.get_atoms()])
            coords2_b = np.array([atom.coord for atom in res2_b.residue.get_atoms()])

            min_len = min(len(coords1_a), len(coords2_a))
            coords1_a = coords1_a[:min_len]
            coords2_a = coords2_a[:min_len]

            min_len = min(len(coords1_b), len(coords2_b))
            coords1_b = coords1_b[:min_len]
            coords2_b = coords2_b[:min_len]

            rmsd_sum += Complex.kabsch_rmsd(coords1_a, coords2_a)
            rmsd_sum += Complex.kabsch_rmsd(coords1_b, coords2_b)

        rmsd_value = rmsd_sum / (len(common_keys) * 2)  # 2 for pairs
        return rmsd_value

# Example usage:
complex1 = Complex('targets/1acb_complex_H.pdb')
complex2 = Complex('/home/as4643/palmer_scratch/Decoys/Supersampled_structures/sampled_1acb/1acb_relaxed/complex.0_13_40_corrected_H_0001.pdb')

print("Interaction Dictionary for Complex 1:")
for key, value in complex1.get_interaction_dict().items():
    print(f"Interaction: {key} - Residue 1: {value[0].get_name()}, Residue 2: {value[1].get_name()}")

print("\nInteraction Dictionary for Complex 2:")
for key, value in complex2.get_interaction_dict().items():
    print(f"Interaction: {key} - Residue 1: {value[0].get_name()}, Residue 2: {value[1].get_name()}")

# Plot interactions for Complex 1
complex1.plot_interactions()

# Calculate RMSD between the two interaction dictionaries
try:
    rmsd = Complex.calculate_rmsd(complex1.get_interaction_dict(), complex2.get_interaction_dict())
    print(f"RMSD between Complex 1 and Complex 2: {rmsd}")
except ValueError as e:
    print(e)
