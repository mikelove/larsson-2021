#!/bin/bash
#
#SBATCH --job-name=download
#SBATCH --time=1440
#SBATCH --mem=1000
#SBATCH --mail-user=milove@email.unc.edu
#SBATCH --mail-type=ALL

wget -o logfile -i fastq_files -nc --random-wait
