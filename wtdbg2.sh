#SBATCH -o wtdbg2%j.out
#SBATCH -e wtdbg2%j.err
#SBATCH -p 2tb
# SBATCH -N 4
# SBATCH -n 64
# SBATCH --cpus-per-task=32
#SBATCH -D /lustre/groups/cbi/ndata
#SBATCH -J wtdbg2
#SBATCH --export=NONE
#SBATCH -t 3-00:00:00
# SBATCH --mem=1k
#SBATCH --mail-type=ALL
#SBATCH --mail-user=azamyatin@email.gwu.edu
# SBATCH --array=1-16
# SBATCH --nice=100

cd /lustre/groups/cbi/ndata

/groups/cbi/wtdbg2/wtdbg2 -x ont \
       -g 280m \
       -t 48 \
       -i ../2018_mopti_inbred_5runs.fq \
       -fo coluzzii
