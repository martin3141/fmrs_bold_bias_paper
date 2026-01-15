#!/bin/bash

#SBATCH --mail-type END
#SBATCH --ntasks 33 
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1
#
#              d-hh:mm:ss
#SBATCH --time 0-06:00:00
#
#SBATCH --account wilsonmp-mrs-analysis
#SBATCH --qos bbdefault

# mail-type options : NONE, BEGIN, END, FAIL
# qos options       : bbshort, bbdefault, castles

set -e

module purge; module load bluebear

module load bear-apps/2024a
module load R/4.5.0-gfbf-2024a
module load R-bundle-CRAN/2025.06-foss-2024a

echo "Modules loaded starting Rscript"

echo "Simulating data..."
Rscript 01_bednarik_simulation.R

echo "Preprocessing data..."
Rscript 02_preprocess_lb.R

echo "Starting analysis..."
Rscript 03_run_fitting.R

echo "Rscript finished :)"
