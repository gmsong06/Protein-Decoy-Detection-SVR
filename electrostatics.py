"""
from pdb_reader import Protein
Protein1 = Protein("file")
Protein1.get_...()

"""
















"""
parser_1 = argparse.ArgumentParser()
parser_1.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser_1.parse_args()

def electrostatics_calc(pdb_file):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    

def process_pdb_file(args):
    path, filename = args
    if filename.endswith('.pdb') and ("NoH" not in filename):
        pdb_path = os.path.join(path, filename)
        return filename[:-4], distance_checker(pdb_path)
    return None

def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    pool = Pool(cpu_count())
    tasks = [(path, filename) for path in paths for filename in os.listdir(path)]
    results = pool.map(process_pdb_file, tasks)
    pool.close()
    pool.join()

    output_csv = f'{pdb_id}_all_contacts.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'all_contacts'])
        for result in results:
            if result:
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
"""