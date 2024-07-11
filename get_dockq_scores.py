import pandas as pd
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
import Bio.PDB
import os
import sys

with open('missing_pdb_files.txt', 'r') as file:
    lines = file.readlines()

lines = [line.strip() for line in lines]

def get_chain_ids(structure):
    """
    Returns two chain ids for the dimer structure
    """
    atoms = 0
    chain_id_1 = ''
    chain_id_2 = ''
    for model in structure:
        for chain in model:
            if not chain_id_1:
                chain_id_1 = chain.id
            elif chain.id != chain_id_1:
                chain_id_2 = chain.id
                return chain_id_1, chain_id_2
            for residue in chain:
                for atom in chain:
                    atoms +=1 
    return [chain_id_1, chain_id_2]
    
structure2 = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', "/home/ms4688/palmer_scratch/Protein-Decoy-Detection-SVR/capri_targets/Target37.pdb")
structure = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', "/home/ms4688/palmer_scratch/Capri_SuperSampled/sampled_T37/sampled_T37_relaxed/complex.0_0_46_corrected_H_0001.pdb")
print(get_chain_ids(structure))
print(get_chain_ids(structure2))

"""
def main():
    scores = {}
    for file in lines:
        #print(file[-3:])
        #print(os.path.dirname(file))
        #folder_name = os.path.basename(os.path.dirname(f"{file[:-4]}.pdb"))
        #print(folder_name)
        #sys.exit()
        pdb_id = file[-3:]
        file_name = file[:-4]
        type = "random"
        if file_name.startswith('complex'):
            if file_name[9] == '_':
                type = "relaxed"

            if type == "relaxed":
                model_path = f"/home/ms4688/palmer_scratch/Capri_SuperSampled/sampled_{pdb_id}/sampled_{pdb_id}_relaxed/{file_name}.pdb"
            else:
                model_path = f"/home/ms4688/palmer_scratch/Capri_SuperSampled/sampled_{pdb_id}/random_negatives/random_{pdb_id}_relaxed/{file_name}.pdb"
        elif file_name.startswith("random"):
            print(file_name)
            model_path = f"/home/ms4688/palmer_scratch/Capri_SuperSampled/sampled_{pdb_id}/random_negatives/random_{pdb_id}_relaxed/{file_name}.pdb"
        
        native = load_PDB(f"/home/ms4688/palmer_scratch/Protein-Decoy-Detection-SVR/capri_targets/Target{pdb_id[-2:]}.pdb")
        print(pdb_id[-2:])
        model = load_PDB(model_path)

        structure = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', model_path)
    
        chains = get_chain_ids(structure)

        chain_map = {chains[0]: chains[0], chains[1]: chains[1]}
        print(chain_map)
        print(run_on_all_native_interfaces(model, native, chain_map=chain_map)[0])
        dockq_score = round(run_on_all_native_interfaces(model, native, chain_map=chain_map)[0][(chains[0], chains[1])]['DockQ'], 3)
        scores[file_name + f"_{pdb_id}"] = dockq_score
        print(dockq_score)

    # Convert the scores dictionary to a DataFrame
    scores_df = pd.DataFrame(list(scores.items()), columns=['Decoy', 'DockQ'])
    scores_df.to_csv("missing_dockq_scores.csv", index=False)
    print("CSV file has been created successfully.")

if __name__ == "__main__":
    main()

"""