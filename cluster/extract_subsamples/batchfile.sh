#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=45g
#SBATCH --job-name=prep_hdf5.49.2
#SBATCH --export=NONE
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=1-5 #%0-12658%200
#unset SLURM_EXPORT_ENV
#export OMP_NUM_THREADS=1

module load gcc/6.2.0
module load python/3.7.4
#module load gsl/2.3 openblas/0.2.19
source /n/groups/reich/hringbauer/explore_ntbk/jptvenv37/bin/activate
module load samtools
module load bcftools

path_vcf="/n/data1/hms/genetics/reich/1000Genomes/ali/WholeGenomeImputation/imputed_r2/v51.1/chr${SLURM_ARRAY_TASK_ID}.bcf"
path_out="/n/groups/reich/hringbauer/git/hapBLOCK/notebook/vignette/data/example_hazelton_chr${SLURM_ARRAY_TASK_ID}.vcf"

bcftools view -s I12439,I12440,I12438,I12896,I21390,I30300 $path_vcf > $path_out

