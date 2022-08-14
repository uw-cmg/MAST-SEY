#!/bin/bash
#
#SBATCH --partition=mem_0096
#SBATCH --qos=mem_0096
#SBATCH -J ptfe-spa
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=400

cd ..
python3 do_all.py
tail -n 1 -q out/*.txt | sort -g > out/sey.plot
date
