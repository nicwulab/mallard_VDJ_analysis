from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.import_pose import CIF_file as CIF
from pyrosetta import init, pose_from_file
import pandas as pd
import os
from Bio import SeqIO
import pandas as pd 


# 2. Initialize the Rosetta scoring & fragment libraries
def pose_from_file(filename):
    """
    Load a pose from a PDB file.
    """
    init(extra_options='-mute all')
    pose = pose_from_file(filename)
    # Run DSSP to get secondary structure
    dssp = Dssp(pose)
    secstruct = dssp.get_dssp_reduced_IG_as_L_secstruct()
    # get the chain infor
    pdb_info = pose.pdb_info()
    Chain = "".join([pdb_info.chain(i+1) for i in range(len(pose.sequence()))])
    TB = pd.DataFrame([list(secstruct), list(Chain), list(pose.sequence())], index= ['SecStruct', 'Chain', 'Sequence']).T
    TB.to_csv(f'result/SecStructure/{filename.split("/")[-1].replace(".cif", ".csv")}', index=False)


# pose_from_file('pdbModel/3jbw.cif')

Flist = [f"AF3/{i}" for i in  os.listdir('AF3/') if "fold_" in i]
List = [[f"{dir}/{i}" for i in os.listdir(dir) if i.endswith('.pdb')] for dir in Flist]

CifList = ['pdbModel/3jbw.pdb']
for i in List:
    CifList.extend(i)

for filename in CifList:
    init(extra_options='-mute all')
    pose = pose_from_file(filename)
    # Run DSSP to get secondary structure
    dssp = Dssp(pose)
    secstruct = dssp.get_dssp_reduced_IG_as_L_secstruct()
    # get the chain infor
    pdb_info = pose.pdb_info()
    Chain = "".join([pdb_info.chain(i+1) for i in range(len(pose.sequence()))])
    TB = pd.DataFrame([list(secstruct), list(Chain), list(pose.sequence())], index= ['SecStruct', 'Chain', 'Sequence']).T
    TB.to_csv(f'result/SecStructure/{filename.split("/")[-1].replace(".pdb", ".csv")}', index=False)


# Read the Aligned Structure


seq_file = 'result/TMalign/fold_rag1_rag2_chicken_model_1.pdb.fasta'

TMalign = 'result/TMalign/'
SeqList = [i for i in os.listdir(TMalign) if i.endswith('.fasta')]

tbList = []
for seq_id in SeqList:
    seq_file = os.path.join(TMalign, seq_id)
    seq_id = seq_file.split('/')[-1].split('.')[0].replace('fold_', '').replace('_model_', '_')
    seqs = []
    for seq_record in SeqIO.parse(seq_file, "fasta"):
        seqs += [str(seq_record.seq)]
    tmp_set = pd.DataFrame([list(i)[:len(seqs[0])] for i in seqs], index = ['human', seq_id]).T
    # Add sequence index
    seq_idx = []
    N = 0
    for i in tmp_set[seq_id].values:
        if i != '-':
            N += 1
            seq_idx.append(N)
        else:
            seq_idx.append(-1)
    # Add sequence index to the dataframe
    tmp_set[seq_id+"_id"] = seq_idx
    tmp_set = tmp_set[~tmp_set.human.isna()]
    tmp_set = tmp_set[tmp_set.human != '-']
    tmp_set['resi'] = tmp_set[seq_id] + tmp_set[seq_id+"_id"].astype(str)
    # Read structure file
    ID = seq_file.split('/')[-1].split('.')[0]
    secFile = f'result/SecStructure/{ID}.csv'
    tb_secStr = pd.read_csv(secFile)
    tb_secStr['resi'] = tb_secStr.Sequence + (tb_secStr.index +1).astype(str)
    tb_secStr = tb_secStr[tb_secStr.Chain == 'A']
    tmp_set.resi.isin(tmp_set.resi)
    tmp_result = pd.merge(tmp_set, tb_secStr[['resi', 'SecStruct']], on='resi', how='left')
    # filter
    tmp_result[seq_id] = tmp_result.resi
    tmp_result[seq_id+"_sec"] = tmp_result.SecStruct
    tmp_result = tmp_result[['human', seq_id,seq_id+"_sec"]]
    tbList.append(tmp_result)

tb_Result = pd.concat([tbList[0]['human']] + [i.iloc[:,1:] for i in tbList], axis = 1)
tb_Result.to_csv('result/Aligned_SecStruct.csv', index=False)









