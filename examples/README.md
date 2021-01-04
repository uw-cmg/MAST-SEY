We will start with preparing the files neccessary for simulation: cumulative integrals of elastic and inelastic cross sections.

Preparing the files from a large energy and q-dependent ELF (file `elf.in`) is a computationally intensive and time consuming process.
Lets divide the preparation to 6 chunks. We start with copying the initial directory into 6 folders.

```bash
for i in {1..6}
do
    cp -r Cu_DFT Cu_DFT-$i
done
```

Then we prepare the input file for parallel execution. The energy (`-e`) ranges from 10 to 1000 eV on a logarithmic scale. The integration (`-i`) of inelastic differential cross sections is carried on 1000 divisions on the energy and 500 divisions on the momentum. Sumrules (`-sumr`) will be calculated to check the sanity of ELF. The elastic cross sections are calculated with the use of Fermi distribution of nuclear charge, Dirac-Fock (DF) model for electron distribution and Furness-McCarthy (FM) exchange potential. A user-provided ELF (calculated with DFT in this case) will be provided, therefore a `-qdep CUSTOM` flag will be set.

```bash
echo "cd Cu_DFT-1; mast_sey prepare -e 10.0 21.2011702065 40 -i 1000 500 -sumr -elastic F DF FM -qdep CUSTOM > Cu_DFT-1.stdout.txt
cd Cu_DFT-2; mast_sey prepare -e 21.6136459698 45.8234586988 40 -i 1000 500 -sumr -elastic F DF FM -qdep CUSTOM > Cu_DFT-2.stdout.txt
cd Cu_DFT-3; mast_sey prepare -e 46.7149692107 99.0412013428 40 -i 1000 500 -sumr -elastic F DF FM -qdep CUSTOM > Cu_DFT-3.stdout.txt
cd Cu_DFT-4; mast_sey prepare -e 100.968080601 214.064146224 40 -i 1000 500 -sumr -elastic F DF FM -qdep CUSTOM > Cu_DFT-4.stdout.txt
cd Cu_DFT-5; mast_sey prepare -e 218.228834836 462.670667132 40 -i 1000 500 -sumr -elastic F DF FM -qdep CUSTOM > Cu_DFT-5.stdout.txt
cd Cu_DFT-6; mast_sey prepare -e 471.672077655 1000.0 40 -i 1000 500 -sumr -elastic F DF FM -qdep CUSTOM > Cu_DFT-6.stdout.txt" > para_prepare
```

if you have 6 cores available you can just execute this file with xargs (replace `-P 0` with the number of cores):

```bash
cat para_prepare | xargs -I EXE -P 0 bash -c EXE
```

or you can use GNU Parallel to distribute the jobs for you:

```bash
parallel < para_prepare
```

you can monitor the progress by looking into the `Cu_DFT-\*.stdout.txt` files.

After the calculations are finished we can gather the results in a new directory, ready for the MC simulation step:

```bash
mkdir Cu_DFT-MC
cat Cu_DFT-[1-6]/elastic\*.in > Cu_DFT-MC/elastic.in
cat Cu_DFT-[1-6]/inelastic\*.in > Cu_DFT-MC/inelastic.in
cat Cu_DFT-[1-6]/energies\*.in > Cu_DFT-MC/energies.in
```

We will also copy the `elf.in`, `material.in` and `dos.in` files.

```bash
cp elf.in Cu_DFT-MC
cp material.in Cu_DFT-MC
cp dos.in Cu_DFT-MC
```

The `mfp\*` files can be also combined to conveniently plot mean free paths.

```bash
cat Cu_DFT-[1-6]/mfp\* > Cu_DFT-MC/mfp.plot
```

We are now ready to perform the MC simulation step, lets move to the newly created directory:

```bash
cd Cu_DFT-MC
```

and prepare the input file for parallel execution. We will run the simulation for energies (`-e`) from 25 eV to 975 eV in 50 eV increments, with 100000 iterations (`-m`) for each energy. Secondaries will be generated from the joint dos (`-dos`). First execution, for initial energy of 25 eV will be run separately, this way a `jdos.in` file will be prepared to be used in subsequent calculations:

```bash
mast_sey -e 25 -m 100000 -dos > 25-Cu_DFT.stdout.txt
```

then, we will prepare an input tile for parallel execution of the MC simulations

```bash
for e in {50..975..50}
do
    echo "mast_sey -e ${e} -m 100000 -dos > ${e}-Cu_DFT.stdout.txt"
done > para_mc
```

and execute that file, either with xargs (replace `-P 0` with your number of cores):

```bash
cat para_prepare | xargs -I EXE -P 0 bash -c EXE &
```

or GNU Parallel:

```bash
parallel < para_prepare &
```
