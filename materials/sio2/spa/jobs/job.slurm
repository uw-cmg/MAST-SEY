#!/bin/bash
#
#SBATCH --partition=mem_0384
#SBATCH --qos=mem_0384
#SBATCH -J sio2-spa
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=400

python3 do_all.py
tail -n 1 -q *txt | sort -g > sey.plot
date
