#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36      # Cores per node
#SBATCH --partition=MT          # Partition Name
##
#SBATCH --job-name=glass
#SBATCH -o test.%N.%j.out         # STDOUT
#SBATCH -e test.%N.%j.err         # STDERR
##

hostname
date

mpirun -np $SLURM_NTASKS lmp_new < lammps_script >& out
