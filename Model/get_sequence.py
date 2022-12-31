"""
Contains a class that takes a single PDB file for the scoring script. It imports the sequence and extracts the C3' atoms in a list type.
The parameter file_path is supposed to be set to the PDB-file of the RNA to be scored.
"""

COL_NAMES = ["Type", "SerNum", "Atom", "ResName", "ChainID", "ResSeqNum", "X", "Y", "Z"]

class Sequence:

    def __init__(self, file_path):
        self.file_path = file_path
        self.seq = []

    def pdb_import(self) ->list:

        with open(self.file_path, "r") as pdbfile:
            for line in pdbfile:
                if line[:4] == "ATOM":       
                    splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21:22], line[22:26], line[30:38], line[38:46], line[46:54]]
                    col= [l.strip() for l in splitted_line]
                    self.seq.append(col)

        if not self.seq:
            raise AttributeError("PDB file could not be imported..\nCheck file format!")
        else:
            print("PDB file imported successfully.")

    def only_c3(self) -> list:
        onlyc3 = [l for l in self.seq if l[2] == "C3\'"]
        self.seq = onlyc3
        print("Imported sequence of ", len(self.seq), " C3\'")