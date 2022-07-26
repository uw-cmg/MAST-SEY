
<img src="misc/MAST-SEY_logo_sm.png" width="50%">


MAST-SEY is an open-source Monte Carlo code capable of predicting secondary electron emission using input data generated entirely from first principle (density functional theory) calculations. It utilises the complex dielectric function and Penn's theory for inelastic scattering processes, and relativistic Schrödinger theory by means of partial-wave expansion method to govern elastic scattering. It allows to not only use the momentum independent (q=0) dielectric function but also to include explicitly calculated momentum dependence, as well as to utilise first-principle density of states in secondary electron generation.  

For more details please refer to the paper which this code accompanies, which is to be cited whenever the code is used:

Maciej P. Polak and Dane Morgan, *MAST-SEY: MAterial Simulation Toolkit for Secondary Electron Yield. A Monte Carlo Approach to Secondary Electron Emission Based On Complex Dielectric Functions*, Comput. Mater. Sci. **193** (2021), 110281 (https://doi.org/10.1016/j.commatsci.2021.110281)

## Installation

1. Download the source code and compile. Optimal performance is achieved with Intel compilers:
```bash
icc -std=c++11 -g -O3 -o mast_sey mast_sey.cpp
```
Compiling with `gcc` is also possible, although the code works twice slower, version > 6.2 of `gcc` is required:
```bash
g++ -std=c++11 -g -O3 -o mast_sey mast_sey.cpp
```
`-Ofast` and architecture specific optimization flags may be used but the performance is unlikely to get noticeably better.

2. [Download the elsepa code](https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/w4hm5vymym-1.zip) (https://doi.org/10.1016/j.cpc.2004.09.006, https://doi.org/10.1016/j.cpc.2020.107704).

3. Unzip the downloaded `w4hm5vymym-1.zip` file and unpack the file inside:
```bash
unzip w4hm5vymym-1.zip
unzip elsepa-2020.zip
```
4. Move the `elsepa2020.patch` patch inside the `elsepa-2020` directory and apply it:
```bash
cd elsepa-2020
patch < ../elsepa2020.patch
```
If for some reason you wish to use the old version of `elsepa` the `elsepa_old.patch` is available too.

5. Compile the patched `elsepa`:
```bash
ifort -o elsepa elscata.f
```
or
```bash
gfortran -o elsepa elscata.f
```
6. Add write permissions to all density files for convenience:
```bash
chmod +w database/z_*
```
7. Make sure that the necessary files are executable:
```bash
chmod +x /path/to/your/elsepa/elsepa /path/to/your/mastsey/mast_sey /path/to/your/mastsey/getDDCS
```
8. Add the directories containing compiled elsepa and mast_sey to your PATH.
```bash
export PATH=/path/to/your/elsepa:/path/to/your/mastsey:${PATH}
```    

## Usage

The code is executed in two steps:
### The "prepare" step
This step postprocesses the input files to a form convenient for the second step to use. It takes the dielectric function `eps.in` or the energy loss function `elf.in`, and using the parameters contained in `material.in`, prepares the cumulative integrals of cross sections. These results are stored in `inelastic.in` and `elastic.in`. Additionally a file `mfp.plot` is generated, and allows for a convenient plotting of the inelastic and elastic mean free paths, which are generated in this step as well. This step is performed only once for each case.

Example `material.in` file for metals containing (the atomic number, volume (A^3), the Fermi energy (eV), the work function (eV)):
```
79
16.929
9.0
5.1
```
Example `material.in` file for insulators containing (the atomic number, volume (A^3), the band gap energy (eV), the width of the valence band (eV), the electron affinity (eV)):
```
79
16.929
6.7
15.8
1.3
```

The `elf.in` file with the custom q-dependency has the following structure:
```
energy[1] ELF[1][q[1]] 
energy[2] ELF[2][q[1]] 
... ...
energy[n] ELF[n][q[1]] 
q[1] q[1]
99999999 99999999
energy[1] ELF[1][q[2]] 
energy[2] ELF[2][q[2]] 
... ...
energy[n] ELF[n][q[2]] 
q[2] q[2]
99999999 99999999
.....
energy[1] ELF[1][q[k]] 
energy[2] ELF[2][q[k]] 
... ...
energy[n] ELF[n][q[k]] 
q[k] q[k]
```

Example `ph.in` file for insulators containing ([energy loss (eV)], eps_0, eps_inf) to include the electron-phonon interaction:
```
0.062 0.15
3.84
2.25
```

The command below is an example of how to run the "prepare" step:
```bash
mast_sey prepare -e 1000 100 -i 100 50 -qdep SPA -elastic P DHFS FM
mast_sey prepare -e 1000 100 -i 100 50 -qdep SPA -elastic P DHFS FM -ins
mast_sey prepare -e 1000 100 -i 100 50 -qdep SPA -elastic P DHFS FM -ins -ph
```
The user should be greeted with the default MAST-SEY output screen. It contains the basic info along with a short feedback on the chosen options and files used. If a basic error is detected, it will be displayed here. In all the input values are correct, a progress bar on the bottom should start filling up (although for accurate calculations it may take a while for even the first bar to appear).

A detailed list of options is given upon execution of the code with the `-h` flag:
```
MAterials Simulation Toolkit for Secondary Electron Emission (MAST-SEY)
Cite as: https://doi.org/10.1016/j.commatsci.2020.XXXXXX
(c) 2020 Maciej P. Polak (mppolak@wisc.edu) & Dane Morgan

Input options:

"prepare" as first argument will run input preparation from "eps/elf.in" and "material.in"
otherwise, the "simulate" version will be executed

"prepare" options:
-ins     calculate properties for insulators
-e       [iniE(eV,optional) range(eV) grid] energy range and grid (def: 1000 500)
-lin     generate the energy grid on a linear scale (default is logarithmic)
-i       [ICS q-int] grids for ICS and q integration (def: 1000 100)
-qdep    [SPA/SSPA/CUSTOM] specify type of q-dependence of ELF (def: SPA)
-sumr    output sum rules for plotting
-saveq   [E_grid q_grid q_max] save q-dependence for plotting
-elastic [nuclear electron exchange (SOLID LDA opt.)] models to use in elastic scattering (def: F TFM FM)
         nuclear: [P]oint/[U]niform/[F]ermi
         electron: [TFM]Thomas–Fermi–Moliere/[TFD]Thomas-Fermi-Dirac/[DHFS]Dirac–Hartree–Fock–Slater/[DF]Dirac-Fock
         exchange: [NO]/[FM]Furness–McCarthy/[TF]Thomas-Fermi/[RT]Riley–Truhlar
         (optional): [SOLID] muffin-tin model potential
         (optional): [LDA] LDA correlation–polarization potential model

"simulate" options:
-ins     simulate for insulators
-e       [incident_energy(eV)] energy of incident energy
-m       [number_of_e-] number of incident electrons (def: 1000)
-core    [energy(eV)] allow secondaries to come from bound states
-dos     [FEG (optional)] generate secondaries from joint DOS from prepared "jdos.in"
         or from parabolic free electron gas approximation (FEG)
-pa      [angle(deg)] angle of incident electrons with respect to surface normal
-coord   save travel paths of e-
-distr   save distribution of secondaries
-noang   use classical approach to inelastic angle scattering
-noout   supress all output

-v       display version of the code
-h       this message


Please be careful when giving input arguments, there is no extensive input checks
Example executions:
mast_sey prepare -e 1000 100 -i 100 50 -elastic F DHFS FM -qdep SPA
mast_sey -e 350 -m 10000

mast_sey prepare -e 1000 100 -i 100 50 -elastic F DHFS FM -qdep SPA -ins
mast_sey -e 350 -m 10000 -ins

mast_sey prepare -e 1000 100 -i 100 50 -elastic F DHFS FM -qdep SPA -ins -ph
mast_sey -e 350 -m 10000 -ins -ph
```

The [examples directory](https://github.com/uw-cmg/MAST-SEY/tree/master/examples) contains examples that showcase most of the capabilities of the code.
They serve as a tutorial for the code and also allow to reproduce the results presented in the paper which this code accompanies:
Maciej P. Polak and Dane Morgan, *MAST-SEY: MAterial Simulation Toolkit for Secondary Electron Yield. A Monte Carlo Approach to Secondary Electron Emission Based On Complex Dielectric Functions*, Comput. Mater. Sci. **193** (2021), 110281 (https://doi.org/10.1016/j.commatsci.2021.110281)

Two material systems are considered: [copper] and [aluminium].

Each directory contains all the neccessary input files and a set of commands that should be executed. It also contains a directory with the calculations already completed, for comparison.

It is highly recommended to execute the code on multiple cores for efficiency. For that users might find particularily helpful the [GNU Parralel](https://www.gnu.org/software/parallel/) software.

