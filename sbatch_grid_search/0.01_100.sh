#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=hp_tuning
##SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-core=1
#SBATCH -c 2
##SBATCH -Q
#SBATCH --mem-per-cpu=50G
##SBATCH --gpus=1
##SBATCH --constraint="p100|v100|rtx2080ti|rtx5000"
#SBATCH -t 48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gmsong1978@gmail.com
#SBATCH --output=hp_tuning.out
#SBATCH --requeue
##SBATCH --array=1-2

# This command navigates to the directory where you submitted the job
cd $SLURM_SUBMIT_DIR
module load miniconda
conda activate SVR

# This command sends the tasklist to all in the array
python /home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/capri_svr_manual_tune.py --remove none --output_path /home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/predictions_capri_grid_search/0.01_100.csv --gamma 0.01 --c 100 --fig_output_name 0.01_100
