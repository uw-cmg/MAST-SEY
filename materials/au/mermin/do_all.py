import os
import numpy as np

energies = [50, 100, 200, 500, 1000, 2000, 4500]
tey = []
sey = []

for energy in energies:
    print(energy)
    os.system("mast_sey -e {} -m 10000 -core 339 356 549 646 763 > ttt".format(energy))
    res = np.loadtxt('ttt', comments='#')
    tey.append(res[1])
    sey.append(res[2])

data = np.column_stack([energies, tey, sey])
np.savetxt('results', data, fmt=['%d','%.2f','%.2f'])
