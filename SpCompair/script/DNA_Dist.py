#!/usr/bin/env python3

pdbFile = 'AF3/fold_rag1_rag2_human/fold_rag1_rag2_human_model_3.pdb'

# Read PDB file by Biopython and calculate the pairwise residues distance between the first two chains.
# Using the SCIPY library for PDB parsing and distance calculation.
from Bio.SeqUtils import seq1
from Bio.PDB import PDBParser
import numpy as np

pdbFile = 'AF3/fold_rag1_rag2_human/fold_rag1_rag2_human_model_3.pdb'

Dirs = [i for i in os.listdir("AF3") if 'fold' in i]

PDBlist = []
for Dir in Dirs:
    PDBlist += [f"AF3/{Dir}/{i}" for i in os.listdir("AF3/" + Dir) if i.endswith('.pdb')]


# Load the PDB file
for pdb_file in PDBlist:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    # Extract atoms from chain A and chain B
    resis_chain_A = []
    resis_chain_B = []
    for model in structure:
        for chain in model:
            if chain.id == "A":
                resis_chain_A.extend([residue for residue in chain])
            elif chain.id == "H":
                resis_chain_B.extend([residue for residue in chain])
    # Convert to numpy arrays of coordinates
    Dist_all = []
    for resiA in resis_chain_A:
        Dist_resi = []
        coords_A = np.array([atom.coord for atom in resiA])
        for resiB in resis_chain_B:
            coords_B = np.array([atom.coord for atom in resiB])
            dis_n = np.min(np.linalg.norm(coords_A[:, np.newaxis, :] - coords_B[np.newaxis, :, :], axis=-1))
            Dist_resi.append(dis_n)
        Dist_all.append(Dist_resi)
    TB = pd.DataFrame(Dist_all, columns = [resiB.get_resname()[1:] + str(resiB.id[1]) for resiB in resis_chain_B], index = [seq1(resiA.get_resname()) + str(resiA.id[1]) for resiA in resis_chain_A])
    TB.to_csv("result/ChainDist/" + pdb_file.split('/')[-1].replace('.pdb', '_distG.csv'), index=True, header=True)
#  translate 3 codon residue name into 1 latter
