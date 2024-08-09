import os
import gzip
import shutil
import argparse


def unzip_pdb_files(input_folder, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file_name in os.listdir(input_folder):
        if file_name.endswith('.pdb.gz'):
            input_file_path = os.path.join(input_folder, file_name)
            output_file_path = os.path.join(output_folder, file_name[:-3])

            with gzip.open(input_file_path, 'rb') as f_in:
                with open(output_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            print(f'Unzipped {file_name} to {output_file_path}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unzip .pdb.gz files")
    parser.add_argument("input_folder", help="Folder path containing .pdb.gz files")
    parser.add_argument("output_folder", help="Folder path to save unzipped .pdb files")

    args = parser.parse_args()

    unzip_pdb_files(args.input_folder, args.output_folder)
