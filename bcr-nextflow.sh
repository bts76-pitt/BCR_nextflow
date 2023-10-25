#!/bin/bash
#SBATCH --job-name=<bcr-nextflow>
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-12:00 #runtime in D-HH:MM
#SBATCH --cpus-per-task=4 # Request that ncpus be allocated per process.
#SBATCH --mem=64
#SBATCH --output=./slurm_log/bcr-nf_%A.out
#SBATCH -e ./slurm_log/bcr-nf_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=your_email@sample.com

module load singularity/3.9.6 nextflow/23.04.2 squashfs-tools/4.4 gcc/12.2.0 r/4.3.0

cd /ix/drajasundaram/drajasundaram/shared_bts76_dhr11/BCR_pipeline/bcr-test
nextflow run ../BCR_nextflow/main.nf -profile slurm -resume -work-dir ./work
