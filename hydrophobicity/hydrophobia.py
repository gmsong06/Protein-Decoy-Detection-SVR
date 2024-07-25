from pdb_reader import Protein
import numpy as np
import pandas as pd
import os
import argparse
import csv
from collections import deque, defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()


def get_residue_name(protein, residue_id):
    for model in protein.structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == residue_id:
                    return residue.resname
    return None

def external_contacts(protein):
    coordinates = []
    atom_names = []
    chains = []
    residue_ids = []

    for model in protein.structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chains.append(chain.id)
                    atom_names.append(atom.get_name())
                    coordinates.append(atom.coord)
                    residue_ids.append(residue.id[1])

    def hydrogen_list():
        for i in range(len(atom_names)):
            first_letter = atom_names[i][0]
            if first_letter == "H":
                atom_names[i] = "hydrogen"
            elif first_letter.isnumeric():
                if atom_names[i][1] == "H":
                    atom_names[i] = "hydrogen"
                else:
                    atom_names[i] = "not_hydrogen"
            else:
                atom_names[i] = "not_hydrogen"

    hydrogen_list()

    dist_thresh = 5  # Distance threshold for contacts
    
    listA = []
    listB = []
    residue_ids_A = []
    residue_ids_B = []
    lst = protein.get_chain_ids()
    for i in range(len(coordinates)):
        if ((chains[i] == lst[0]) and (atom_names[i] == "not_hydrogen")):
            listA.append(coordinates[i])
            residue_ids_A.append(residue_ids[i])
        elif ((chains[i] == lst[1]) and (atom_names[i] == "not_hydrogen")):
            listB.append(coordinates[i])
            residue_ids_B.append(residue_ids[i])      

    res_in_contact = []
    
    for i,posA in enumerate(listA):
        for j,posB in enumerate(listB):
            if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                cont = [residue_ids_A[i], residue_ids_B[j]]
                if cont not in res_in_contact:
                    res_in_contact.append(cont)
    
    return listA, listB, residue_ids_A, residue_ids_B, res_in_contact

def get_hydro_hits(protein):

    hydrophobicity_dict = {
            "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
            "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
            "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
            "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
        }

    hitsA = []
    hitsB = []

    listA, listB, resA, resB, res_list = external_contacts(protein)

    # print(f"Contacts in Chain A: {len(resA)} + Contacts in Chain B: {len(resB)}")
    # if len(resA) == len(resB):
    #     print("Lengths are equal. Proceeding to compare hydrophobicities.")
    # print("-----------------------------------------------------")

    for i in range(len(resA)):
        resnameA = get_residue_name(protein, resA[i])
    for i in range(len(resB)):
        resnameB = get_residue_name(protein, resB[i])
        
        if abs(hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]) <= 0.2:
            hitsA.append(resA[i])
            hitsB.append(resB[i])


    return hitsA, hitsB

def internal_contacts(protein):
    listA, listB, resA, resB, res_list = external_contacts(protein)
    
    lstA = {}
    lstB = {}
    
    for i,posA in enumerate(listA):
        tempA=[]
        for j,posB in enumerate(listA):
            if i != j:
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= 3.5:
                    tempA.append(resA[j])
        lstA[i] = tempA
    
    for i,posA in enumerate(listB):
        tempB=[]
        for j,posB in enumerate(listB):
            if i != j:
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= 3.5:
                    tempB.append(resB[j])
        lstB[i] = tempB

    return resA, resB, lstA, lstB
            
def create_graph(protein):
    hydroA = {}
    hydroB = {}
    
    resA, resB, graphA, graphB = internal_contacts(protein)
    listA, listB = get_hydro_hits(protein)

    for i in range(len(resA)):
        tempA =[]
        if resA[i] in listA:
            tempA.append(True)
        hydroA[i] = tempA

    for i in range(len(resB)):
        tempB =[]
        if resB[i] in listB:
            tempB.append(True)
        hydroB[i] = tempB

    return (graphA, graphB, hydroA, hydroB)
    # print(graphA)
    # print(graphB)
    # print(hydroA)
    # print(hydroB)

# prot = Protein("/Users/smriti/Desktop/aeop/Protein-Decoy-Detection-SVR/targets/1c3a_complex_H.pdb")
#create_graph(prot)

lst = [(1, [1]), (2, [2, 1]), (3, [3, 2]), (4, [3, 2, 2]), (5, [4, 3, 5, 2])]

def score_fnc(lst):
    max_dist = len(lst)
    avg = 0
    score = 0
    SD = 0
    ns = []
    for dist in lst:
        dist_allowed, islands = dist
        tot = 0
        for island in islands:
            tot += len(islands) * island
        ns.append(tot)
        avg += tot
        score += tot/dist_allowed

    avg = avg/len(lst)
    for val in ns:
        SD += ((val - avg)**2)
    
    SD = SD/(max_dist-1)
    weighted_avg = score/max_dist

    return score, weighted_avg, avg, SD


def bfs(adj_list, start_node, visited, distances):
    q = deque([start_node])
    visited[start_node] = True
    distances[start_node] = 0

    while q:
        current_node = q.popleft()

        for neighbor in adj_list[current_node]:
            if not visited[neighbor]:
                visited[neighbor] = True
                distances[neighbor] = distances[current_node] + 1
                q.append(neighbor)

def compute_distances(adj_list):
    distances = defaultdict(dict)
    
    for start_node in adj_list:
        visited = {node: False for node in adj_list}
        node_distances = {node: float('inf') for node in adj_list}
        bfs(adj_list, start_node, visited, node_distances)
        
        for node, dist in node_distances.items():
            if dist != float('inf'):
                distances[start_node][node] = dist
    
    return distances


def find_islands(adj_list, dist_allowed, hydro_reaction):
    distances = compute_distances(adj_list)
    islands = []

    in_island = set()

    visited = {}

    for node in adj_list:
        visited[node] = False
    
    for node in adj_list:
        # Only begin exploring from nodes that are hydro and aren't already in an island
        if not hydro_reaction[node] or node in in_island:
            continue
        
        # print(f"Starting node is {node}")
        # We start with a hydro node so the distance starts at 0

        current_island = [node]
        in_island.add(node)

        q = deque()
        q.append((node, 0))

        visited = set()

        while q:
            # print(f"Length of q is {len(q)}")
            curr_node_info = q.popleft()

            curr_node = curr_node_info[0]
            dist_from_hydro_node = curr_node_info[1]

            visited.add(curr_node)

            # print(f"Current node is {curr_node}")

            for dist_node in distances[curr_node]:
                dist_from_hydro_node = curr_node_info[1]
                # print(f"Dist node is {dist_node}")
                # print(f"Dist from hydro node is {dist_from_hydro_node}")
                distance = distances[curr_node][dist_node]
                # print(f"Node is {dist_node}. Distance is {distance}")

                # Check if it can be explored
                if distance + dist_from_hydro_node <= dist_allowed and dist_node not in visited and dist_node not in in_island:
                    # print(f"{dist_node} node made it into q")
                    

                    # If the new node is a hydro node
                    if hydro_reaction[dist_node]:
                        # print(f"{dist_node} node is a hydro reaction so distance from hydro node is 0")
                        dist_from_hydro_node = 0
                        current_island.append(dist_node)
                        in_island.add(dist_node)
                        # print(f"Here is the current island: {current_island}")
                    else:
                        # print(f"{dist_node} is not a hydro reaction so the distance from hydro node is {dist_from_hydro_node + distance}")
                        dist_from_hydro_node += distance
                    
                    q.append((dist_node, dist_from_hydro_node))
            # print()

        islands.append(current_island)

    return [len(island) for island in islands]


def get_max_dist(adj_list, hydro_reaction):

    def bfs_dist(adj_list, start_node, hydro_reaction, visited, distances):
        q = deque([start_node])
        visited[start_node] = True
        distances[start_node] = 0

        while q:
            current_node = q.popleft()

            for neighbor in adj_list[current_node]:
                if not visited[neighbor]:
                    visited[neighbor] = True
                    distances[neighbor] = distances[current_node] + 1
                    q.append(neighbor)

    max_distance = 0
    nodes_with_reaction = [node for node in hydro_reaction if hydro_reaction[node]]

    for start_node in nodes_with_reaction:
        # Initialize visited and distance dictionaries
        visited = {node: False for node in adj_list}
        distances = {node: float('inf') for node in adj_list}

        bfs_dist(adj_list, start_node, hydro_reaction, visited, distances)

        # Calculate maximum distance to nodes with hydro_reaction as True
        for end_node in nodes_with_reaction:
            if distances[end_node] != float('inf'):
                max_distance = max(max_distance, distances[end_node])

    return max_distance


def get_final_island_data(adj_list, hydro_reaction):
    max_dist = get_max_dist(adj_list, hydro_reaction)

    final_data = []

    for i in range(max_dist):
        dist_allowed = i + 1

        islands = find_islands(adj_list, dist_allowed, hydro_reaction)

        final_data.append((dist_allowed, islands))

    return final_data

def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    for path in paths:
        print(f"Path is {path}")
        for filename in os.listdir(path):
            print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                print(f"Processing {filename}")
                prot = Protein(pdb_path)
                results.append((filename[:-4], (get_residues(prot))))
            else:
                print(f"File did not pass requirements.")

def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-4:]
            print(f"PDB id is {pdb_id}")
            process_pdb_folder(full_folder_path, pdb_id)
            print("DONE----------------------------------------------------------------------")

    graph = {
        0: [1, 4, 6],
        1: [0, 2, 4, 7],
        2: [1, 3, 8],
        3: [2, 5, 9],
        4: [0, 1],
        5: [3, 6],
        6: [0, 5],
        7: [1, 8],
        8: [2, 7, 9],
        9: [3, 8]
    }


    hydro_reaction = {
        0: True, 
        1: False, 
        2: False,
        3: False,
        4: True,
        5: True,
        6: True,
        7: False,
        8: True,
        9: False,
    }
    
    dist_allowed = 2
    find_islands(graph, dist_allowed, hydro_reaction)

    print(get_final_island_data(graph, hydro_reaction))
    print(score_fnc(get_final_island_data(graph, hydro_reaction)))
    # print(islands)
    # print(compute_distances(graph))

if __name__ == "__main__":
    main(args.pdb_folder)
