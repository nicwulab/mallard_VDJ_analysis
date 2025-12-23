#!/usr/bin/env python3


'''
cd-hit -i final/IGHV.fa -o cd-hit/IGHV.fa -s 0.75
cd-hit -i final/IGLV.fa -o cd-hit/IGLV.fa -s 0.75
'''


import os
from Bio import SeqIO
from Bio.Seq import Seq



def IGV_Clu(CD_hit, input_fa, output_fa):
    with open(CD_hit, 'r') as f:
        clusters = f.read().split('>Cluster ')[1:]

    cls = clusters[0]  # Example: first cluster
    ID_uniq = []
    ID_group = {}
    for cls in clusters:
        IDs = [i.replace('>', '').replace('.', '') for i in cls.split() if '>' in i]
        ID_1st = IDs[0]
        ID_uniq.append(ID_1st)
        ID_group[ID_1st] = IDs[1:]

    ID_value = {}
    for ID in ID_uniq:
        ID_value[ID] = int(ID.split('-')[-1])

    # sort by value
    ID_sorted = sorted(ID_value.items(), key=lambda x: x[1], reverse=False)

    ID_rename = {}
    N = 0
    for ID,_ in ID_sorted:
        N += 1
        new_ID =  ID.split('-')[0] + '-' + str(N)
        ID_rename[ID] = new_ID
        N_sup = 1
        if len(ID_group[ID]) > 0:
            ID_rename[ID] = f"{new_ID}*01"
            for ID_g in ID_group[ID]:
                N_sup += 1
                ID_rename[ID_g] = f"{new_ID}*{str(N_sup).zfill(2)}"

    with open(input_fa, 'r') as f_in, open(output_fa, 'w') as f_out:
        Old_file = f_in.read()
        for key, value in ID_rename.items():
            Old_file = Old_file.replace(f">{key} ", f">{value} ")
        f_out.write(Old_file)


if __name__ == '__main__':
    # IGHV
    CD_hit = 'cd-hit/IGHV.fa.clstr'
    input_fa = 'final/IGHV.fa'
    output_fa = 'final/IGHV_cls.fa'
    IGV_Clu(CD_hit, input_fa, output_fa)
    # IGLV
    CD_hit = 'cd-hit/IGLV.fa.clstr'
    input_fa = 'final/IGLV.fa'
    output_fa = 'final/IGLV_cls.fa'
    IGV_Clu(CD_hit, input_fa, output_fa)
