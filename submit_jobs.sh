#!/bin/bash

# Define the main directory containing the sbatch files and result folders
main_dir="/vast/palmer/scratch/ohern/sr2562"
# Define the main folder containing subfolders to be run
main_folder="${main_dir}/Supersampled-Structures/groups"

# Iterate over each subfolder in the main folder
for subfolder in "${main_folder}"/*/; do
    # Get the subfolder name
    subfolder_name=$(basename "$subfolder")
    
    # Update the sbatch script with the new paths, job name, and output file
    sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=A${subfolder_name}|" "$sbatch_hydro"
    sed -i "s|^#SBATCH --output=.*|#SBATCH --output=${main_dir}/${subfolder_name}/A${subfolder_name}.out|" "$sbatch_hydro"
    
    # Update any custom path inside the sbatch script (if applicable)
    # Assuming there is a line in the sbatch script that specifies the working directory
    # For example, "#SBATCH --chdir=/path/to/directory"
    sed -i "s|^#SBATCH --chdir=.*|#SBATCH --chdir=${subfolder}|" "$sbatch_file"
    
    # Submit the modified sbatch script
    sbatch "$sbatch_hydro"

done
