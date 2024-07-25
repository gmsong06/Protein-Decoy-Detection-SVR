#!/bin/bash
#SBATCH --partition=pi_ohern,gpu
#SBATCH --job-name=T46
##SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=1
#SBATCH -c 2
##SBATCH -Q
#SBATCH --mem-per-cpu=50G
##SBATCH --gpus=1
##SBATCH --constraint="p100|v100|rtx2080ti|rtx5000"
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gmsong1978@gmail.com
#SBATCH --output=T46.out
#SBATCH --requeue
##SBATCH --array=1-2

# This command navigates to the directory where you submitted the job
cd $SLURM_SUBMIT_DIR
module load miniconda
conda activate SVR

# This command sends the tasklist to all in the array
python -m hydrophobicity.hydro_capri /home/as4643/palmer_scratch/Decoys/Capri_SuperSampled/Capri_subsets/T46
