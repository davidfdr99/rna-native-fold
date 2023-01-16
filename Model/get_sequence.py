"""
Contains functions to preprocess single PDB file for scoring. It imports the sequence and extracts the C3' atoms in a list type.
The parameter file_path is supposed to be the path of a PDB-file of the RNA to be scored.
"""

COL_NAMES = ["Type", "SerNum", "Atom", "ResName", "ChainID", "ResSeqNum", "X", "Y", "Z"]

def pdb_import(file_path: str) ->list:
    """
    Imports PDB file and 
    """

    seq = []

    with open(file_path, "r") as pdbfile:
        for line in pdbfile:
            if line[:4] == "ATOM":       
                splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21:22], line[22:26], line[30:38], line[38:46], line[46:54]]
                col= [l.strip() for l in splitted_line]
                seq.append(col)

    if not seq:
        raise AttributeError("PDB file could not be imported..\nCheck file format!")

    print("PDB file imported successfully.")
    return seq

def only_c3(seq) -> list:
    onlyc3 = [l for l in seq if l[2] == "C3\'"]
    print("Imported sequence of ", len(onlyc3), " C3\'")
    return onlyc3

def get_sequence(file: str) -> list:
    """
    Wrapper that calls pdb_import and only_c3 subsequently
    """

    seq = pdb_import(file_path=file)
    c3_seq = only_c3(seq)

    return c3_seq

