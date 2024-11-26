#!/bin/bash

#SBATCH -J job
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -o %j.out
#SBATCH -e %j.err

module purge

# Shell 脚本任务提交
bash script.sh

# Python 脚本任务提交
conda activate <virtual_env>
# source path/activate <virtual_env>

python script.py
