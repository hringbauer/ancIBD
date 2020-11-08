#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8g
#SBATCH --job-name=simulate_ibd
#SBATCH --export=NONE
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=0-4%5  #%0-12658%200
#unset SLURM_EXPORT_ENV
#export OMP_NUM_THREADS=1

module load gcc/6.2.0
module load python/3.7.4
#module load gsl/2.3 openblas/0.2.19
source /n/groups/reich/hringbauer/explore_ntbk/jptvenv37/bin/activate

# Execute the following tasks
python3 run_single.py $SLURM_ARRAY_TASK_ID 
