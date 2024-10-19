#!/bin/bash

#SBATCH -J vasp_local
#SBATCH -p 64c512g
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -o %j.out
#SBATCH -e %j.err

module purge

module load intel-oneapi-compilers/2021.4.0
module load intel-oneapi-mpi/2021.4.0
module load intel-oneapi-mkl/2021.4.0

ulimit -s unlimited
ulimit -l unlimited

mpirun ${HOME}/yangsl/bin/vasp.5.std
