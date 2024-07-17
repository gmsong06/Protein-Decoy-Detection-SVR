#!/bin/bash

# Define the main directory containing the sbatch files and result folders
main_dir="/vast/palmer/scratch/ohern/sr2562"
# Define the main folder containing subfolders to be run
main_folder="${main_dir}/Supersampled_structures/groups"

# Iterate over each subfolder in the main folder
for subfolder in "${main_folder}"/*/; do
    # Get the subfolder name
    subfolder_name=$(basename "$subfolder")
    
    # Update the sbatch script with the new paths, job name, and output file
    sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=C${subfolder_name}|" "$sbatch_hydro"
    sed -i "s|^#SBATCH --output=.*|#SBATCH --output=${main_dir}/${subfolder_name}/C${subfolder_name}.out|" "$sbatch_hydro"
    
    # Update the line where the python script is run
    # Replace "my_script.py" with your actual script name and "argument" with your argument pattern
    sed -i "s|^python /vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/hydrophobe.py .*|python /vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/hydrophobe.py ${main_dir}/${subfolder_name}|" "$sbatch_hydro"


    
    # Submit the modified sbatch script
    sbatch "$sbatch_hydro"

done
