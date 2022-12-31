import os
import argparse

from Model.training_model import get_freq_matrix, get_observed_frequency, get_reference_frequency, get_pseudo_energy

FILE_DIR = './FilesPDB'

def make_pseudoenergy_dir():
    """
    Creates a directory named 'Pseudoenergy' where the 10 files for each base pair combination
    containing the pseudoenergy for each distance bin in the interval [0, 20] Å will be stored.
    """

    current = os.path.dirname(__file__)
    dir = os.path.join(current, 'Pseudoenergy')
    
    if not os.path.isdir(dir):
        os.makedirs(dir, mode= 0o777)
        print(f"Successfully created directory: \n{dir}\n")
    else:
        print(f"Pseudoenergy files will be stored in existing folder: \n{dir}.\n")

def get_energy_files(pseudo_energy: dict):
    """
    Writes a file containing the pseudoenergy as for each base pair combination.
    The files have 20 lines corresponding to the pseudo energy for each distance bin (in the interval [0, 20] Å)
    """

    make_pseudoenergy_dir() # Creates the Pseudoenergy directory where the files will be stored in

    for key, values in pseudo_energy.items():
        with open(f"./Pseudoenergy/{key}.txt", 'w') as profile:
            for value in values:
                profile.write(f"{value}\n")
         
def main(arg):
    """
    Creates the model that can thereafter be used for RNA PDB-files scoring.
    Argument is the directory of the PDB-Files with which the model will be trained. The default is './FilesPDB'.
    The functions are imported from the 'training model' directory.
    """
    
    ff, ss = get_freq_matrix(dir_path=arg)
    total = 0
    for value in ss.values():
        total += value

    # Get observed frequencies
    obs_freq = get_observed_frequency(ff, ss)
    
    # Get reference frequency
    ref_freq = get_reference_frequency(ff, total)

    # Get pseudo-energy
    u = get_pseudo_energy(obs_freq, ref_freq)

    # Write a file for each pseudo-energy profile:
    get_energy_files(u) 

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument('files', type=str, nargs='?', default=FILE_DIR, help="You can specify a folder containing a PDB-files library for model training.")
    args=parser.parse_args()

    main(arg=args.files)
