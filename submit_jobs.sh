#!/bin/bash

# Define the main directory containing the sbatch files and result folders
main_dir="/vast/palmer/scratch/ohern/sr2562"
# Define the main folder containing subfolders to be run
main_folder="${main_dir}/Supersampled_structures/groups"

# Define the sbatch file
sbatch_hydro="${main_dir}/sbatch_hydro.sh"

# Check if the sbatch file exists
if [[ ! -f "$sbatch_hydro" ]]; then
    echo "Error: sbatch file not found at $sbatch_hydro"
    exit 1
fi

# Iterate over each subfolder in the main folder
for subfolder in "${main_folder}"/*/; do
    # Get the subfolder name
    subfolder_name=$(basename "$subfolder")

    echo "Processing subfolder: $subfolder_name"

    # Update the sbatch script with the new paths, job name, and output file
    sed -i "s|^#SBATCH --job-name=.*|#SBATCH --job-name=C${subfolder_name}|" "$sbatch_hydro"
    sed -i "s|^#SBATCH --output=.*|#SBATCH --output=${main_dir}/${subfolder_name}/C${subfolder_name}.out|" "$sbatch_hydro"
    sed -i "s|^#SBATCH --error=.*|#SBATCH --error=${main_dir}/${subfolder_name}/C${subfolder_name}.err|" "$sbatch_hydro"

    # Update the line where the python script is run
    sed -i "s|^python /vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/hydrophobe.py .*|python /vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/hydrophobe.py ${subfolder}|" "$sbatch_hydro"

    # Verify changes
    echo "Updated sbatch file for $subfolder_name:"
    grep "^#SBATCH --job-name=" "$sbatch_hydro"
    grep "^#SBATCH --output=" "$sbatch_hydro"
    grep "^#SBATCH --error=" "$sbatch_hydro"
    grep "^python /vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/hydrophobe.py" "$sbatch_hydro"

    # Submit the modified sbatch script
    sbatch "$sbatch_hydro"
    echo "Submitted sbatch file for $subfolder_name"
done
