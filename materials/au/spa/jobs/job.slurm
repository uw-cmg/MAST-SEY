#!/bin/bash
#
#SBATCH --partition=mem_0384
#SBATCH --qos=mem_0384
#SBATCH -J au-spa
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=30

cd ..
python3 do_all.py
tail -n 1 -q out-feg/*txt | sort -g > out-feg/sey.plot
date
