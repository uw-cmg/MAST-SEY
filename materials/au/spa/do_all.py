import os
from multiprocessing import Pool

energies = [50, 100, 200, 500, 1000, 2000, 4500]

def run_mast(i):
	print(energies[i])
	os.system("mast_sey -e {} -m 1000 -core 549 > {}-spa.stdout.txt".format(energies[i],energies[i]))


if __name__ == '__main__':
	with Pool(processes=4) as pool:
		pool.map(run_mast, range(len(energies)))
