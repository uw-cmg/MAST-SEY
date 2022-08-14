#!/bin/bash
#
#SBATCH --partition=mem_0384
#SBATCH --qos=mem_0384
#SBATCH -J mast_prep
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>

mast_sey prepare -e 4900 4900 -i 1000 500 -qdep SPA -elastic F DHFS FM
