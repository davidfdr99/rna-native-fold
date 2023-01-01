import os
import math
import argparse
from datetime import datetime

from Model.get_sequence import Sequence
from Model.training_model import get_distances
from Model.testing_model import to_output_dir, write_results_file, uniquify, linear_interpolation

RNA_TEST = './Test/4gxy.pdb'
TRAINING_DIR = './Pseudoenergy'
OUTPUT_DIR='./Output'

def main(test, training, output):
    """
    First, the distances between residues are computed with the same restrictions as in model training. The result will be printed to the terminal
    and simultaneously stored in an output text file to save it.
    """

    # Import the PDB-file
    sample = Sequence(test)
    sample.pdb_import()
    sample.only_c3()

    # Create distances dictionary for the sample
    dist_sample = get_distances(sample.seq)

    # Compute Gibb's free energy
    gibbs_energy = linear_interpolation(sample_distances=dist_sample, training_dir=training)

    #Write output text file
    to_output_dir(output)
    sample_name = test.strip('.pdb').split('/')[-1]
    file_path = output + "/rnfold_" + sample_name + ".tsv"
    if os.path.isfile(file_path):
        file_path = uniquify(file_path) # Adds counter to path so files are not overwritten
    with open(file_path, 'w') as f:
        write_results_file(f, gibbs_energy, sample_name) 
    
    print(f"The estimated Gibb's free energy is {gibbs_energy}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('test_rna', type=str, nargs='?', default=RNA_TEST, help='Please parse a valid path to the PDB file of the RNA sequence to be tested!')
    parser.add_argument('training', type=str, nargs='?',default=TRAINING_DIR, help='Please parse a valid path to the directory containing pseudoenergy profile files computed by the \`training.py\` script!')
    parser.add_argument('output', type=str, nargs='?', default=OUTPUT_DIR, help="You can choose a output directory where the result will be stored. Please parse a vaild path. If the folder does not exist yet, the directory will be created.")
    args=parser.parse_args()

    main(test=args.test_rna, training=args.training, output=args.output)