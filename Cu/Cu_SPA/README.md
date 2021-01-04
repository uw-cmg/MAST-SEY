We will start with preparing the files neccessary for simulation: cumulative integrals of elastic and inelastic cross sections.

The energy (`-e`) ranges from the Fermi energy to 1000 eV on a logarithmic scale and distributed over 250 points. The integration (`-i`) of inelastic differential cross sections is carried on 1000 divisions on the energy and 500 divisions on the momentum. Sumrules (`-sumr`) will be calculated to check the sanity of ELF. The elastic cross sections are calculated with the use of Fermi distribution of nuclear charge, Dirac-Fock (DF) model for electron distribution and Furness-McCarthy (FM) exchange potential. A Simple Penn Approximation will be used, so the `-qdep SPA` flag will be set.

```bash
mast_sey prepare -e 1000 250 -i 1000 500 -sumr -elastic F DF FM -qdep SPA > Cu_SPA-prep.stdout.txt &
```

this step should take no more than 5 minutes. You can monitor the progress by looking into the `Cu_SPA-prep.stdout.txt` files.

We are now ready to perform the MC simulation step, lets prepare the input file for parallel execution. We will run the simulation for energies (`-e`) from 25 eV to 975 eV in 50 eV increments, with 100000 iterations (`-m`) for each energy. Secondaries will be generated from the joint dos (`-dos`). First execution, for initial energy of 25 eV will be run separately, this way a `jdos.in` file will be prepared to be used in subsequent calculations:

```bash
mast_sey -e 25 -m 100000 -dos > 25-Cu_SPA.stdout.txt &
```

This should only take around a minute, similarly as before, you can monitor the progress in the output txt file. Then, we will prepare an input file for parallel execution of the MC simulations:

```bash
for e in {975..75..-50}
do
    echo "mast_sey -e ${e} -m 100000 -dos > ${e}-Cu_SPA.stdout.txt"
done > para_mc
```

and execute that file, either with xargs (replace `-P 0` with your number of cores):

```bash
cat para_mc | xargs -I EXE -P 0 bash -c EXE &
```

or GNU Parallel:

```bash
parallel < para_prepare &
```

this step may take up to an hour. For practice, you may reduce the number of iterations by an order of magnitude. As before, you can monitor the progress in the output txt files.

Once all the calculations are done, the results can be gathered from all the outut files and the SEY curve can be plotted with any plotting software:

```bash
tail -n 1 -q *txt | sort -g > sey.plot
```
