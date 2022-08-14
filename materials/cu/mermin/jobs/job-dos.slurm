#!/bin/bash
#
#SBATCH --partition=mem_0096
#SBATCH --qos=mem_0096
#SBATCH -J cu-m-dos
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=40

cd ..
python3 do_all_dos.py
tail -n 1 -q out-dos/*txt | sort -g > out-dos/sey.plot
date
