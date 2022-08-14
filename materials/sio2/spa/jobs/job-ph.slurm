#!/bin/bash
#
#SBATCH --partition=mem_0384
#SBATCH --qos=mem_0384
#SBATCH -J sio2-spa-ph
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=240

cd ..
python3 do_all_ph.py
tail -n 1 -q out-ph/*.txt | sort -g > sey-ph.plot
date
