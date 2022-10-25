#!/bin/bash
#
#SBATCH --job-name=snake
#SBATCH --time=1440
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

snakemake -j 15 --latency-wait 30 --cluster 'sbatch -n 12 --mem=10000 --time=10'
