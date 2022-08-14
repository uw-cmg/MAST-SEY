#!/bin/bash
#
#SBATCH --partition=mem_0096
#SBATCH --qos=mem_0096
#SBATCH -J cu-m-u
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=40

cd ..
python3 do_all.py
tail -n 1 -q out-u/*txt | sort -g > out-u/sey.plot
date
