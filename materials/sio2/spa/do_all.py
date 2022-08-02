import os
from multiprocessing import Pool

energies = [50, 100, 200, 300, 500, 700, 1000, 2000, 3000, 4500]

def run_mast(i):
	print(energies[i])
	os.system("mast_sey -e {} -m 10000 -ins > {}-spa.stdout.txt".format(energies[i],energies[i]))


if __name__ == '__main__':
	with Pool(processes=10) as pool:
		pool.map(run_mast, range(len(energies)))

