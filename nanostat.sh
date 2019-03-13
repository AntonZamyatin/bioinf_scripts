#!/bin/bash

#SBATCH -o NanoStat%j.out
#SBATCH -e NanoStat%j.err
#SBATCH -p short
#SBATCH -N 1
# SBATCH -n 16 
# SBATCH --cpus-per-task=16
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH  -J NanoStat
#SBATCH --export=NONE
#SBATCH -t 2-00:00:00
# SBATCH --mem=0
#SBATCH --mail-user=azamyatin
# SBATCH --array=1-16
#SBATCH --nice=100

module load python3/3.5.2
export PY_USER_BIN=$(python -c 'import site; print(site.USER_BASE + "/bin")')
export PATH=$PY_USER_BIN:$PATH
NanoStat --fastq 2018_mopti_inbred_5runs.fq -t 16 -o out -n stat
