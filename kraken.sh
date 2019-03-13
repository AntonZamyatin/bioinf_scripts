#!/bin/bash

#SBATCH -o kraken%j.out
#SBATCH -e kraken%j.err
#SBATCH -p short
#SBATCH -N 1
# SBATCH -n 16
# SBATCH --cpus-per-task=16
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH  -J kraken
#SBATCH --export=NONE
#SBATCH -t 1:00:00
# SBATCH --mem=0
#SBATCH --mail-type=ALL
#SBATCH --mail-user=azamyatin@email.gwu.edu
# SBATCH --array=1-16
# SBATCH --nice=100

module load gcc/8.1.0
cd /lustre/groups/cbi/ndata
kraken/kraken2 --db KrakenDB_ALL assembled/miniasm/miniasm.fa \
               --output miniasm_out > reports/kraken_out
