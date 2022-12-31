An Objective Function for the RNA folding problem
============
This programme calculates the native fold for a given RNA chain, being the one with the lowest Gibb's free energy. 

This is computed by interpolating interatomic distances over a given dataset of distance distributions trained from experimentally determined (thus known) 3D structures.

The project was developed as coursework for the *Structural Bioinformatics of RNA* course of the Master 2 GENIOMHE of *Université Paris Saclay - Université d'Evry Val d'Essonne* under the supervision of Prof. Guillaume Postic @guipostic.

---
## Usage 

To get the source code, clone this repository:
```bash
$ git clone https://github.com/davidfdr99/rna-native-fold.git
```

 `script3.py` can be run in a terminal. It takes as argument the path of a PDB file with the sample RNA sequence:

```bash
$ python script3.py FilesPDB/4gxy.pdb
```
The default sequence is Adenosylcobalamin riboswitch (PDB: [DOI 10.2210/pdb4GXY/pdb](https://10.2210/pdb4GXY/pdb)).

---
## Model training

The folder `ModelTraining` contains the scripts with all functions to compute the estimated Gibb's free energy. 

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

For each combination pair, `script1.py` generates a file with one score per distance bin. `script2.py` creates corresponding plots. 

---
## Interaction Profiles

The `InteractionProfiles` folder contains the estimated Pseudo-Energy profiles for each base pair combination. The maximum score is set to 10. 

---
## License [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The source code of this programme is licensed under the MIT license.