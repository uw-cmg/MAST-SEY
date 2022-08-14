#!/bin/bash
#
#SBATCH --partition=mem_0384
#SBATCH --qos=mem_0384
#SBATCH -J pmma-spa-prep
#SBATCH -N 1
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=<olgaridzel@yandex.ru>
#SBATCH --time=01:00:00

cd ..
mast_sey prepare -e 4900 3000 -i 500 500 -qdep SPA -elastic F DHFS FM -ins -ph
