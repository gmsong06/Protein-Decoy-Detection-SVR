#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "<gamma> <c> <output_file> <fig_log_name> <submit>"
    exit 1
fi

# Assign the argument to a variable
GAMMA=$1
C=$2
OUTPUT_FILE=$3
FIG_LOG_NAME=$4
SUBMIT=$5

FULL_OUTPUT_FILE="/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/predictions_capri_grid_search/${OUTPUT_FILE}.csv"

TEMPLATE="/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/sbatch_grid_search/sbatch_grid_search.sh"

NEW_NAME="/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/sbatch_grid_search/${OUTPUT_FILE}.sh"

cp "$TEMPLATE" "$NEW_NAME"

if [ $? -eq 0 ]; then
    echo "File copied to: $NEW_NAME"
else
    echo "Failed to copy file"
    exit 1
fi

echo "python /home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/capri_svr_manual_tune.py --remove none --output_file ${FULL_OUTPUT_FILE} --gamma ${GAMMA} --c ${C} --fig_output_name ${FIG_LOG_NAME}" >> "${NEW_NAME}"

if [ $? -eq 0 ]; then
    echo "Line added to file: $NEW_NAME"
else
    echo "Failed to add line to file"
    exit 1
fi

if [ "$SUBMIT" == "y" ]; then
    # Execute sbatch command on the copied file
    cd "/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR"

    sbatch "$NEW_NAME"

    # Check if the sbatch command was executed successfully
    if [ $? -eq 0 ]; then
        echo "sbatch command executed successfully on $NEW_NAME"
    else
        echo "Failed to execute sbatch command"
        exit 1
    fi
else
    echo "Did not submit."
fi

exit 0
