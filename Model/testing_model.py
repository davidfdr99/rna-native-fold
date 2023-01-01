"""
This file contains several functions used in the rnascore.py script
"""

import os
import math
from datetime import datetime

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