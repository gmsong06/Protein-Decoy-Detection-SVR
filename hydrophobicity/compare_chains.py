from pdb_reader import Protein
import hydrophobicity.hydrophobia as hydrophobia
import numpy as np
import pandas as pd
import os
import argparse
import csv
import math

# parser = argparse.ArgumentParser()
# parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
# args = parser.parse_args()

def score_fnc(islandsA, islandsB, hydro_hitsA, hydro_hitsB):
    # create dictionaries to map residues to their respective islands
    islandsA_dict = {res: i for i, island in enumerate(islandsA) for res in island}
    islandsB_dict = {res: i for i, island in enumerate(islandsB) for res in island}
    
    total_matches = 0
    total_residues = 0
    patch_scores = []
    
    for islandA in islandsA:
        island_match_count = 0
        
        for resA in islandA:
            if resA in hydro_hitsA:
                index = hydro_hitsA.index(resA)
                resB = hydro_hitsB[index]
                
                if resB in islandsB_dict:
                    islandB_idx = islandsB_dict[resB]
                    
                    # check other residues in the same islandA
                    matching_residues = [hydro_hitsB[hydro_hitsA.index(other_resA)] 
                                         for other_resA in islandA if other_resA in hydro_hitsA]
                    
                    # count matches in the same islandB
                    match_count = sum(1 for r in matching_residues if islandsB_dict.get(r, -1) == islandB_idx)
                    
                    if match_count > 0:
                        island_match_count += match_count
        
        island_size = len(islandA)
        patch_score = island_match_count / island_size if island_size > 0 else 0
        patch_scores.append(patch_score)
        total_matches += island_match_count
        total_residues += island_size

    # calculate the percentage of residues in a patch on chain a that match with residues in the same patch on chain b
    patch_alignment_score = total_matches / total_residues if total_residues > 0 else 0
    
    # calculate the overall percentage of patches on chain a that have consistent matches with patches on chain b
    consistent_patches = sum(1 for score in patch_scores if score > 0.5)
    chain_patch_consistency_score = consistent_patches / len(islandsA) if len(islandsA) > 0 else 0
    tot = 0
    for score in patch_scores:
        tot += score
    avg_patch_score = tot/len(patch_scores)
    return patch_alignment_score, chain_patch_consistency_score, avg_patch_score

prot = Protein("/Users/smriti/Desktop/aeop/Protein-Decoy-Detection-SVR/complex.0_0_11_corrected_H_0001.pdb")
graphA, graphB, hydroA, hydroB = hydrophobia.create_graph(prot)
islandsA, final_island_dataA = hydrophobia.get_final_island_data(graphA, hydroA)
islandsB, final_island_dataB = hydrophobia.get_final_island_data(graphB, hydroB)
hydro_hitsA, hydro_hitsB = hydrophobia.get_hydro_hits(prot)
score = score_fnc(islandsA, islandsB, hydro_hitsA, hydro_hitsB)
print(score)

'''
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
                graphA, graphB, hydroA, hydroB = hydrophobia.create_graph(prot)
                islandsA, final_island_dataA = hydrophobia.get_final_island_data(graphA, hydroA)
                islandsB, final_island_dataB = hydrophobia.get_final_island_data(graphB, hydroB)
                hydro_hitsA, hydro_hitsB = hydrophobia.get_hydro_hits()
                score = score_fnc(islandsA, islandsB, hydro_hitsA, hydro_hitsB)

                results.append((filename[:-4], score))
            else:
                print(f"File did not pass requirements.")

    output_csv = f'/home/as4643/palmer_scratch/hydro_results/islands/{pdb_id}_hydro_islands.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(["pdb_file", "patch_alignment_score", "chain_patch_consistency_score", "avg_patch_score"])
        for result in results:
            writer.writerow(result)


def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-4:]
            print(f"PDB id is {pdb_id}")
            process_pdb_folder(full_folder_path, pdb_id)
            print("DONE----------------------------------------------------------------------")

    

if __name__ == "__main__":
    main(args.pdb_folder)
'''