#!/bin/bash

#SBATCH -J vasp_sy
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o %j.out
#SBATCH -e %j.err

module purge

module load vasp/5.4.4-intel-2021.4.0

ulimit -s unlimited
ulimit -l unlimited

mpirun vasp_std
