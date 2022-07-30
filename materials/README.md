An example script for the parallel MC simulation once the prepare step has been done.

```python
import os
from multiprocessing import Pool

num_cores = 4

energies = [50, 100, 200, 500, 1000, 2000, 4500]

def run_mast(i):
	print(energies[i])
	os.system("mast_sey -e {} -m 1000 > e{}.stdout.txt".format(energies[i],energies[i]))


if __name__ == '__main__':
	with Pool(processes=num_cores) as pool:
		pool.map(run_mast, range(len(energies)))
```
