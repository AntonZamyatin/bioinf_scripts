#!/bin/bash

#SBATCH -o miniasm%j.out
#SBATCH -e miniasm%j.err
#SBATCH -p 2tb
#SBATCH -N 1
# SBATCH -n 16 
# SBATCH --cpus-per-task=16
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH  -J miniasm
#SBATCH --export=NONE
#SBATCH -t 16:00:00
# SBATCH --mem=0
#SBATCH --mail-user=azamyatin@email.gwu.edu
# SBATCH --array=1-16
#SBATCH --nice=100

/groups/cbi/miniasm/miniasm -f 2018_mopti_inbred_5runs.fq reads.paf.gz > miniasm.gfa
