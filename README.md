
<img src="MAST-SEY_logo_sm.png" width="50%">

MAST-SEY is an open-source Monte Carlo code capable of predicting secondary electron emission using input data generated entirely from first principle (density functional theory) calculations. It utilises the complex dielectric function and Penn's theory for inelastic scattering processes, and relativistic Schrödinger theory by means of partial-wave expansion method to govern elastic scattering. It allows to not only use the momentum independent (q=0) dielectric function but also to include explicitly calculated momentum dependence, as well as to utilise first-principle density of states in secondary electron generation. 

## Installation

1. Download the source code and compile with `gcc`, version > 6.2 of `gcc` is recommended.
```bash
g++ -std=c++11 -g -O3 -o mast_sey mast_sey.cpp
```
2. Download the elsepa code: https://md-datasets-cache-zipfiles-prod.s3.eu-west-1.amazonaws.com/5zzrz874tt-1.zip (https://doi.org/10.1016/j.cpc.2004.09.006)
3. Apply the `elscata.patch` in the downloaded directory
```bash
patch elscata.patch elscata.f
```
4. Copy the executables `mast_sey1` and `elscata` to a convenient location
5. Add that location to your PATH:
```bash
export PATH=${PATH}:/complete/path/to/your/mast_sey
```    
    - you can add that line to your .bashrc of you dont want to execute it each time
6. Make sure that the copied files are executable:
```bash
chmod +x mast_sey elsepa
```

## Usage

The code is executed in two steps:
1. "prepare" stage:
This step postprocesses the input files to a form convenient for the second step to use. It takes the dielectric function `eps.in` or the energy loss function `elf.in`, and using the parameters contained in `material.in`, prepares the cumulative integrals of cross sections. These results are stored in `inelastic.in` and `elastic.in`. Additionally a file `mfp.plot` is generated, and allows for a convenient plotting of the inelastic and elastic mean free paths, which are generated in this step as well. This step is performed only once for each case.

The command below is an example of how to run the "prepare" step:
```bash
mast_sey -e 1000 100 -i 100 50 -qdep DFT -elastic P DHFS FM
```
In this example, the cumulative integrals will be calculated in the energy range (`-e`) between the Fermi energy (read from `material.in`) and 1000 eV, on a grid of 100 points on a logarithmic scale. The integrals in energy and q (momentum) space will be calculated on grids of 100 and 50 equally spaced points, respectively. A DFT calculated q-dependence of energy loss function will be used. The inelastic scattering cross ections will be calculated using a model with [P]oint nuclear charge, [D]irac-[H]artree-[F]ock-[S]later electron model, and [F]urness-[M]cCarthy exchange.
The `examples` directory containes all the necceary input files to run this example.
## Contributing


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
