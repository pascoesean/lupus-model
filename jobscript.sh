#!/bin/bash
#SBATCH -A b1042             ## Allocation
#SBATCH -p genomics               ## Queue
#SBATCH --job-name=gctp_stc_hg
#SBATCH -t 12:00:00             ## Walltime/duration of the job
#SBATCH -N 1                    ## Number of Nodes
#SBATCH --mem=40G               ## Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --ntasks-per-node=6     ## Number of Cores (Processors)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=seanpascoe2024@u.northwestern.edu

module purge
eval "$(conda shell.bash hook)"
conda activate stomics-env
python3 gsd.py