#!/bin/bash

inpfile=$1
sjobfile=$(basename $inpfile .inp).sjob

ncpus=$2

cat > $sjobfile <<EOF
#!/bin/bash
#SBATCH --partition=CPU_rune
#SBATCH -J $(basename $inpfile .inp)
#SBATCH -N 1
#SBATCH --ntasks-per-node=$ncpus
#SBATCH --cpus-per-task=1

source /usr/local/_tci_software_environment_.sh
module purge
module load rune/molpro/2022.1

mkdir -p /scratch/\$USER/tmpdir_\$SLURM_JOBID
molpro -n \$SLURM_NTASKS_PER_NODE -d /scratch/\$USER/tmpdir_\$SLURM_JOBID -I \$PWD -W \$PWD -a -t 1 $inpfile
rm -rf /scratch/\$USER/tmpdir_\$SLURM_JOBID
EOF

echo "$sjobfile was created for you! It will use $ncpus CPUs "
