"""
Functions are used for training the model (getting the pseudoenergy over input PDB-files) as well as for scoring.
"""

import math
import os

COL_NAMES = ["Type", "SerNum", "Atom", "ResName", "ChainID", "ResSeqNum", "X", "Y", "Z"]
BASES_DICT = {"AA": [], "AU": [], "AC": [], "AG": [], "UU": [], "UC": [], "UG": [], "CC": [], "CG": [], "GG": []}

def get_interatomic_distance(x1, y1, z1, x2, y2, z2):
    """
    Formula of the distance between two atoms in the sequence.
    """

    d = ((x2 -x1)**2) + ((y2-y1)**2) + ((z2-z1)**2)
    return(math.sqrt(d))

def get_distances(seq: list) -> dict:
    """
    Gets the distance between all atoms in a sequence parsed as list, given the restriction of a pairwise distance
    greater than 3 Å.
    """

    base_comb = {"AA": [], "AU": [], "AC": [], "AG": [], "UU": [], "UC": [], "UG": [], "CC": [], "CG": [], "GG": []}
    
    for i, u in enumerate(seq):
        for j, v in enumerate(seq[i:]):

            if int(v[5]) - int(u[5]) > 3:
                dist = get_interatomic_distance(float(u[6]), float(u[7]), float(u[8]), float(v[6]), float(v[7]), float(v[8]))
                if dist <= 20:
                    comb = u[3] + v[3]
                    if comb in base_comb:
                        base_comb[comb].append(dist)
                    elif comb[::-1] in base_comb:
                        base_comb[comb[::-1]].append(dist)       

    return base_comb

def get_freq_matrix(dir_path: str) ->list:
    """
    Takes as input the directory containing a RNA PDB-file library to train the model. It iterates over all files and then computes the frequency
    matrix for each distance bin. Returns frequency matrix and pairwise frequency sums as dictionaries.
    """

    # First create dictionaries for frequency matrix and combination sums
    freq_matrix = {key: [0 for num in range(21)] for key in BASES_DICT}
    pair_freq_sum = {comb: 0 for comb in BASES_DICT}

    # Read PDB Files one by one:
    print(f'Started importing PDB-files from {dir_path}.')
    for filename in os.listdir(dir_path):

        pdb = []

        f = os.path.join(dir_path, filename)
        # Checking if it is a valid file:
        if os.path.isfile(f) and filename.endswith(".pdb"):
            with open(f, "r") as pdbfile:
                for line in pdbfile:
                    if line[:4] == "ATOM":       
                        splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21:22], line[22:26], line[30:38], line[38:46], line[46:54]]
                        col= [l.strip() for l in splitted_line]
                        pdb.append(col)

                if not pdb:
                    raise AttributeError("PDB file could not be imported..\nCheck file format!")
            
            # Only C3' atoms
            onlyc3 = [l for l in pdb if l[2] == "C3\'"]
            counts = get_distances(onlyc3)

            # Calculate the observed frequency counts and sums for each combination, increase frequency matrix at corresponding distance bin
            for key, values in counts.items():
                if key in freq_matrix:
                    for dist in values:
                        freq_matrix[key][(math.floor(dist))+1] += 1
        
    # Calculate sums of frequency counts for each combination
    for key, values in freq_matrix.items():
        freq_sum = sum(values)
        pair_freq_sum[key] = freq_sum
    
    print("Successfully imported PDB files and caculated frequency matrix")
    
    return freq_matrix, pair_freq_sum

def get_observed_frequency(freq_count: dict, pair_freq_sum: dict) -> dict:
    """
    Computes observed frequencies
    """

    freq_obs = {key: [] for key in freq_count}
    for key, values in freq_count.items():
        for value in values:
            if pair_freq_sum[key] != 0:
                freq_obs[key].append(value/pair_freq_sum[key])
            else:
                freq_obs[key].append(value)
    
    print("Successfully computed observed frequencies")
    
    return freq_obs

def get_reference_frequency(freq_count: dict, total: int) -> dict:
    """
    Computes reference frequency
    """

    ref_freq = [0 for i in range(21)]
    for key in freq_count.keys():
        for i, value in enumerate(freq_count[key]):
            ref_freq[i] += (value/total)

    print("Successfully computed reference frequency")
    
    return ref_freq

def get_pseudo_energy(obs_freq: dict, ref_freq: list) -> dict:
    """
    Computes Pseudoenergy as -log of the ratio of observed to reference frequency.
    """

    pseudo_energy = {key: [] for key in obs_freq}

    for key in obs_freq.keys():
        for i, value in enumerate(obs_freq[key]):
            if value != 0:
                pseudo_energy[key].append(-math.log(value/ref_freq[i]))
            else:
                pseudo_energy[key].append(10)
    
    print("Successfully computed Pseudoenergy profiles")
    
    return pseudo_energy        