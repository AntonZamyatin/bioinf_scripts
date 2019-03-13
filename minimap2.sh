#!/bin/bash

#SBATCH -o minimap%j.out
#SBATCH -e minimap%j.err
#SBATCH -p short
# SBATCH -N 4
# SBATCH -n 64
# SBATCH --cpus-per-task=32
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH -J minimap
#SBATCH --export=NONE
#SBATCH -t 3-00:00:00
# SBATCH --mem=1k
#SBATCH --mail-type=ALL
#SBATCH --mail-user=azamyatin@email.gwu.edu
# SBATCH --array=1-16
# SBATCH --nice=100

cd /groups/cbi/minimap2

minimap2 -ax map-ont reference/chromosomes.fasta \
          2018_mopti_inbred_5runs.fq > all_chromo_alignment.sam
