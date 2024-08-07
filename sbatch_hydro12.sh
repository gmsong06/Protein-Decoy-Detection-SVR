#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=12
##SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=1
#SBATCH -c 2
##SBATCH -Q
#SBATCH --mem-per-cpu=10G
##SBATCH --gpus=1
##SBATCH --constraint="p100|v100|rtx2080ti|rtx5000"
#SBATCH -t 168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gmsong1978@gmail.com
#SBATCH --output=/vast/palmer/scratch/ohern/as4643/hydro_results/patches/out_files/12.out
#SBATCH --requeue
##SBATCH --array=1-2

# This command navigates to the directory where you submitted the job
cd $SLURM_SUBMIT_DIR
module load miniconda
conda init
conda activate SVR

# This command sends the tasklist to all in the array
python -m hydrophobicity.compare_chains /vast/palmer/scratch/ohern/as4643/Decoys/Supersampled_structures/groups/12/
