#!/bin/bash

#SBATCH -o canu%j.out
#SBATCH -e canu%j.err
#SBATCH -p short
#SBATCH -N 10
# SBATCH -n 16
# SBATCH --cpus-per-task=16
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH  -J canu
#SBATCH --export=NONE
#SBATCH -t 24:00:00
# SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=azamyatin@email.gwu.edu
# SBATCH --array=1-16
# SBATCH --nice=100

module load canu
module load jdk/1.8.0
cd /lustre/groups/cbi/ndata
canu -p coluzii \
     -d assembled/canu \
     genomeSize=300m \
     useGrid=false \
     -nanopore-raw 2018_mopti_inbred_5runs.fq > canu_out

