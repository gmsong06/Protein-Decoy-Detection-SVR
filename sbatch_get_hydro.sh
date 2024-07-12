#!/bin/bash
#SBATCH --partition=pi_ohern,gpu
#SBATCH --job-name=get_hydro
##SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=1
#SBATCH -c 2
#SBATCH -Q
#SBATCH --mem-per-cpu=10G
#SBATCH --gpus=1
#SBATCH --constraint="p100|v100|rtx2080ti|rtx5000"
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smritirajan211@gmail.com
#SBATCH --output=get_hydro.out
#SBATCH --requeue
##SBATCH --array=1-2

# This command navigates to the directory where you submitted the job
cd $SLURM_SUBMIT_DIR
module load miniconda
conda activate SVR

# This command sends the tasklist to all in the array - FIX BEFORE RUNNING!!!!!!!!!!
python /home/sr2562/palmer_scratch/Protein-Decoy-Detection-SVR/hydrophobe.py --remove none --output_path /home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/predictions/all.csv
