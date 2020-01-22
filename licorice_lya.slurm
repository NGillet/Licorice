#!/bin/sh

#SBATCH -J Lic_lya 
#SBATCH -t 24:00:00
#SBATCH --mem=118000MB
#SBATCH --cpus-per-task=24
#SBATCH --constraint="HSW24"
#SBATCH --ntasks=500
#SBATCH --nodes=500
# Nombre de noeuds:
#SBATCH --ntasks-per-node=1
#SBATCH --output=licorice.log
#SBATCH --error=licorice.err

module purge
module load intel/17.0
module load intelmpi/2018.1.163 
#module load openmpi/intel/2.0.1
#module load intel/17.2
#module load openmpi/intel/.1.8.8
export FORT_BUFFERED=false

ulimit -s unlimited
ulimit -l unlimited
#ulimit -a unlimited
export OMP_NUM_THREADS=24
export KMP_AFFINITY=granularity=fine,compact,1,0,verbose
export KMP_STACKSIZE=100m
export FORT_FMT_RECL=300
#mpirun -n 32 -perhost 2 ./licorice
srun --mpi=pmi2 -K1 -m block:block -c 24 --resv-ports -n $SLURM_NTASKS ./licorice
#mpirun ./licorice_hybrid.out 
