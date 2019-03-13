#!/bin/bash

#SBATCH -o flye%j.out
#SBATCH -e flye%j.err
#SBATCH -p 256gb
# SBATCH -N 4
# SBATCH -n 64
# SBATCH --cpus-per-task=32
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH -J flye
#SBATCH --export=NONE
#SBATCH -t 3-00:00:00
# SBATCH --mem=1k
#SBATCH --mail-type=ALL
#SBATCH --mail-user=azamyatin@email.gwu.edu
# SBATCH --array=1-16
# SBATCH --nice=100

module load flye
cd /lustre/groups/cbi/ndata
flye --nano-raw 2018_mopti_inbred_5runs.fq \
          --genome-size 300m \
          --threads 16 \
          --out-dir /lustre/groups/cbi/ndata/assembled/Flye > flye_out.txt

