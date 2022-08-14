#!/bin/bash
#
#SBATCH --partition=mem_0384
#SBATCH --qos=mem_0384
#SBATCH -J cu-m-prep
#SBATCH -N 1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>

cd ..
mast_sey prepare -e 4900 3000 -i 500 500 -qdep CUSTOM -elastic F DHFS FM
