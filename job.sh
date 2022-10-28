#!/bin/bash

# Submit the pipeline as a job with srun job.sh

# Modified from https://github.com/mschubert/clustermq/blob/master/inst/LSF.tmpl
# under the Apache 2.0 license:
#SBATCH --job-name=SCISSORS_Sims
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --mem-per-cpu=64G
#SBATCH --cpus-per-task=1

module load r # Comment out if R is not an environment module.
R CMD BATCH run.R

# Removing .RData is recommended.
rm -f .RData
