An Objective Function for the RNA folding problem
============
This programme calculates the native fold for a given RNA chain, being the one with the lowest Gibb's free energy. 

This is computed by interpolating interatomic distances over a given dataset of distance distributions trained from experimentally determined (thus known) 3D structures.

The project was developed as coursework for the *Structural Bioinformatics of RNA* course of the Master 2 GENIOMHE of *Université Paris Saclay - Université d'Evry Val d'Essonne* under the supervision of Prof. Guillaume Postic [@guipostic](https://github.com/guipostic).

---
# Installation 

To get the source code, clone this repository:
```bash
$ git clone https://github.com/davidfdr99/rna-native-fold.git
```
It is recommended to initialise a new virtual environment (e.g. using conda) and then installing all the dependencies:

```bash
$ conda create -n [ENV-NAME]
$ pip install -r requirements.txt
```
---
# Usage

## Fast Start

The default test sequence is Adenosylcobalamin riboswitch (PDB: [DOI 10.2210/pdb4GXY/pdb](https://10.2210/pdb4GXY/pdb)). The PDB file can be found in the `Test` directory. To get the free Gibb's energy according to the model using all default parameters and resulting files, run the `rnascore.py` script in a terminal without any additional arguments.
```bash
$ python3 rnascore.py
```

---
## Running the scripts

To train the model with a custom library of RNA PDB-files, `training.py` generates a file with one score per distance bin for each combination pair. It takes as single terminal argument the path of the directory containing all the PDB files. Default is `./FilesPDB`.
```bash
$ python3 training.py ./FilesPDB
```

`plotting.py` creates corresponding plots. It takes as single argument the directory where the plots are to be saved. Default is `./InteractionProfiles`.
```bash
$ python3 plotting.py ./InteractionProfiles
```

`rnascore.py` is the script used to score a single PDB file. It takes as arguments the path to the PDB file to be tested (test), the path to the training directory containing the 10 files of pseudoenergy profiles (training) as well as a path to an output directory (output).

```bash
$ python3 rnascore.py ./Test/4gxy.pdb ./Pseudoenergy ./Output
```

---
## Model training

The folder `Model` contains the scripts with all functions to compute the estimated Gibb's free energy. 

### Constraints

* Only C3' atoms are taken into account
* Only "intrachain" distances are considered
* Residues must be separated by at least 3 positions on the sequence

### Observed Frequency

The probability (i.e. frequency) of observing two residues *i* and *j* separated by distance bin *r* is calculated as follows:

$$ f_{i,j} ^{OBS}(r) = { N_{i,j}(r) \over N_{i,j} } $$

where N_{i,j(r)} is the count of i and j within the distance bin r, and Ni,j is the count of i and j for
all distance bins. Only distance intervals of 0 to 20 Å are taken into account.

### Reference Frequency

The reference frequency is the same formula except that the different residue types (A, U, C, G) are indistinct ("X"):

$$ f_{X,X} ^{REF}(r) = { N_{X,X}(r) \over N_{X,X}} $$

### Pseudo-Energy

The score for each residue pair is computed as followed:

$$ u_{i,j}(r) = { -log \left( f _{i,j} ^{OBS}(r) \over f_{i,j} ^{REF}(r) \right) } $$

### Linear Interpolation

The final score of the Gibb's free energy for a given RNA is computed by linearly interpolating the pair-wise distances of the priorly computed Pseudoenergy profiles. Therefore, the following formula is applied for each distance and the result appended to the total energy:

$$ y = y_1 + { (x - x_1)(y_2 - y_1) \over x_2-x_1 } $$

Where x is the coordinate to perform interpolation and
y is the interpolated value $^1$.

---
## Options

For the three scripts, a help option is available by calling --help, e.g. 

```bash
$ python3 rnascore.py --help

usage: rnascore.py [-h] [test_rna] [training] [output]

positional arguments:
  test_rna    Please parse a valid path to the PDB file of the RNA sequence to be tested!
  training    Please parse a valid path to the directory containing pseudoenergy profile files computed by the \`training.py\` script!
  output      You can choose an output directory where the result will be stored. Please parse a vaild path. If the folder does not exist yet, the directory will be
              created.

optional arguments:
  -h, --help  show this help message and exit
```
---
## References

$^1$ Linear Interpolation Formula: https://www.geeksforgeeks.org/linear-interpolation-formula/ 

---
## License [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The source code of this programme is licensed under the MIT license.