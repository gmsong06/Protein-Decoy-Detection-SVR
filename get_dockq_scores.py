import pandas as pd
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
import Bio.PDB

with open('missing_pdb_files.txt', 'r') as file:
    lines = file.readlines()

lines = [line.strip() for line in lines]

def get_chain_ids(structure):
    """
    Returns two chain ids for the dimer structure
    """
    chain_id_1 = ''
    chain_id_2 = ''
    for model in structure:
        for chain in model:
            if not chain_id_1:
                chain_id_1 = chain.id
            elif chain.id != chain_id_1:
                chain_id_2 = chain.id
                return chain_id_1, chain_id_2
    return [chain_id_1, chain_id_2]


def main():
    scores = {}
    for file in lines:
        pdb_id = file[-4:]
        file_name = file[:-5]
        type = "random"
        if file_name.startswith('complex'):
            if file_name[9] == '_':
                type = "relaxed"

            if type == "relaxed":
                model_path = f"/home/as4643/palmer_scratch/Decoys/Supersampled_structures/sampled_{pdb_id}/{pdb_id}_relaxed/{file_name}.pdb"
            else:
                model_path = f"/home/as4643/palmer_scratch/Decoys/Supersampled_structures/sampled_{pdb_id}/random_negatives/rand_{pdb_id}_relaxed/{file_name}.pdb"
        elif file_name.startswith("random"):
            print(file_name)
            model_path = f"/home/as4643/palmer_scratch/Decoys/Supersampled_structures/sampled_{pdb_id}/random_negatives/rand_{pdb_id}_relaxed/{file_name[:-len('_corrected_H_0001_')]}.pdb"
        
        native = load_PDB(f"/home/as4643/palmer_scratch/Decoys/Supersampled_structures/targets/{pdb_id}_complex_H.pdb")
        model = load_PDB(model_path)

        structure = Bio.PDB.PDBParser(QUIET=True).get_structure('protein', model_path)
        chains = get_chain_ids(structure)

        chain_map = {chains[0]: chains[0], chains[1]: chains[1]}

        dockq_score = round(run_on_all_native_interfaces(model, native, chain_map=chain_map)[0][(chains[0], chains[1])]['DockQ'], 3)
        scores[file_name + f"_{pdb_id}"] = dockq_score
        print(dockq_score)

    # Convert the scores dictionary to a DataFrame
    scores_df = pd.DataFrame(list(scores.items()), columns=['Decoy', 'DockQ'])
    scores_df.to_csv("missing_dockq_scores.csv", index=False)
    print("CSV file has been created successfully.")

if __name__ == "__main__":
    main()
