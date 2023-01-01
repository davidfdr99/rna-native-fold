import os
import math
import argparse
from datetime import datetime

from Model.get_sequence import Sequence
from Model.training_model import get_distances

RNA_TEST = './Test/4gxy.pdb'
TRAINING_DIR = './Pseudoenergy'
OUTPUT_DIR='./Output'

def to_output_dir(path: str):
    """
    Creates output directory if it does not exist yet.
    """
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            raise("Error whilst creating new output directory. Please choose a valid path or use default option!")
    print(f'Results will be stored in {path}')

def linear_interpolation(sample_distances: dict, training_dir: str):
    """
    Computes the Gibb's free energy of a given RNA sample through linear interpolation with the pseudoenergy profiles of the model.
    It takes as arguments a dictionary containing the distances per nucleotide pair of the sample (sample_distances) as well as the path 
    to the directory containing the pseudoenergy profiles. It return the Gibb's free energy.
    """

    # Creates dictionary per distance bin for each nucleotide pair of the energy profiles
    reference = {}
    for filename in os.listdir(training_dir):
        base_comb = filename.strip('.txt')

        with open(f'{training_dir}/{filename}', 'r') as f:
            distances = []
            for line in f:
                distances.append(float(line.strip('\n')))

            reference[base_comb] = {num: 0 for num in range(21)}
            for i, dist in enumerate(distances):
                reference[base_comb][i] = float(dist)

    #Set Gibb's energy to 0
    gibbs_energy = 0

    #Iterate over sample distances and calculate Gibb's free energy through linear interpolation
    for key, value in sample_distances.items():
            for dist in value:
                
                profile_lower_bound = reference[key][math.floor(dist)]
                profile_upper_bound = reference[key][math.ceil(dist)]

                energy = profile_lower_bound + (dist - math.floor(dist)) * (profile_upper_bound-profile_lower_bound)
                gibbs_energy += energy

    return gibbs_energy

def write_results_file(file, gibbs_energy, sample):
    """
    Function to write output file content
    """

    file.write(f"Sample\t{sample}\n")
    file.write(f"Gibb's Free Energy\t{gibbs_energy}\n")
    file.write(f"Run\t{datetime.date(datetime.now())}\t{datetime.time(datetime.now())}\n")

def uniquify(path):
    filename, extension = os.path.splitext(path)
    counter = 1

    while os.path.exists(path):
        path = filename + " (" + str(counter) + ")" + extension
        counter += 1

    return path

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