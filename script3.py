import os, sys
import math

from Model.get_sequence import Sequence
from Model.training_model import get_distances

RNA_TEST = 'FilesPDB/4gxy.pdb'

# --> dictionary including distance bin for each base_comb in sample

def linear_interpolation(sample_dict):

    gibbs_energy = 0
    training_dir = './InteractionProfiles/Training'

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
    
    for key, value in sample_dict.items():
            for dist in value:

                x1, x2, y1, y2 = math.floor(dist), math.ceil(dist), reference[key][math.floor(dist)], reference[key][math.ceil(dist)]
                energy = y1 + (dist - x1) * (y2-y1)/(x2-x1)

                gibbs_energy += energy

    return gibbs_energy

def main(arg):
    sample = Sequence(arg)
    sample.pdb_import()
    sample.only_c3()

    dist_sample = get_distances(sample.seq)

    gibbs_energy = linear_interpolation(dist_sample)
    
    print(f"The estimated free Gibb's energy is {gibbs_energy}")


if __name__ == "__main__":

    if len(sys.argv) > 1:
        arg = sys.argv[1]
    else:
        arg = RNA_TEST
    main(arg)