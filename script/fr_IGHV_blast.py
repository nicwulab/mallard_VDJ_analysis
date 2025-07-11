#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
This script is used to identify the potential IGHV framework regions from the blast results.
The main pipeline is as follows:
    1. Blast the connected IGHV regions to the genome (fwr1 to fwr3)
    2. read the blast results and remove redundancies

Blast results was done by other scripts, we just need to read the results here.
In the end, we get 78 if the potential IGHV regions. 1 of the sequence is aligned reversely and be removed.

Previous, we use the fwr1 to perform the blast and the we got 75. So, it seems this result is better.
'''

import os
import numpy as np
import pandas as pd
from collections import Counter

fr_list = [i for i in os.popen('ls result/ptg00*.tsv').read().split() if 'fr1_3'  in i][-1:]

fr_all = pd.DataFrame()
for fr_file in fr_list:
    fr = fr_file.split('/')[-1].split('.')[0]
    tb = pd.read_csv(fr_file, sep='\t')
    tb.columns = ['qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps'] 
    tb['seg'] = fr.split('_')[-1]
    fr_all = pd.concat([fr_all, tb])

# split by scalfold and identify the directions
DictV = {}

tb = fr_all.copy()
tb['forward'] = tb.sstart[tb.seg == 'fr1'].mean() - tb.sstart[tb.seg == 'fr4'].mean() < 0
tb['seg_f'] = tb.sstart - tb.send < 0

# remove the reversed 
# tb = tb[tb.forward == tb.seg_f]
# remove duplicates 
tb = tb.sort_values(['sstart', 'send'], ascending=[True, False])
# remove duplicates
tb = tb[~tb.sstart.duplicated()]
tb = tb[~tb.send.duplicated()]
# sort again 
tb = tb.sort_values(['sstart', 'send'], ascending=[True, False])
# remove overlaps
# mask_overlap = np.where(tb.sstart[1:].to_numpy() - tb.send[:-1].to_numpy() < 0)[0]
# tb = tb.iloc[mask, :]
# tb = tb[tb.seg_f == True]

TB_all = tb
TB_all.to_csv('result/PotentialHV.tsv', sep='\t') 

# ## take the seq from the genome
# save the sequences

TB_all.sstart
TB_all.send

## take the seq from the genome
from Bio import SeqIO
seq_loc="../data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

## all qstart = 1, so I don't need to change the start and end
Seq_aa = []
Seq_nr = []
for i in range(TB_all.shape[0]):
    seq_id = f"{TB_all.iloc[i].qacc}:{TB_all.iloc[i].sacc}"
    if TB_all.seg_f.iloc[i]:
        seq_aa = str(seq_dict[TB_all.sacc.iloc[i].split('_')[0]][(TB_all.iloc[i].sstart-1):(TB_all.iloc[i].sstart-1) +390].translate().seq)
        seq_nr = str(seq_dict[TB_all.sacc.iloc[i].split('_')[0]][(TB_all.iloc[i].sstart-1):(TB_all.iloc[i].sstart-1) +390].seq)
    else:
        seq_aa = str(seq_dict[TB_all.sacc.iloc[i].split('_')[0]][(TB_all.iloc[i].send-1):(TB_all.iloc[i].send-1) + 390].translate().seq)
        seq_nr = str(seq_dict[TB_all.sacc.iloc[i].split('_')[0]][(TB_all.iloc[i].send-1):(TB_all.iloc[i].send-1) + 390].seq)
    Seq_aa.append(f">{seq_id}\n{seq_aa}")
    Seq_nr.append(f">{seq_id}\n{seq_nr}")
    if 'CACGGTGGCTCAGAGCCCGTGCGCGGCTGTCACAAAACC' in seq_nr:
        print(i)

with open('result/PotentialHV_aa.fa', 'w') as f:
    f.write('\n'.join(Seq_aa))
with open('result/PotentialHV_nr.fa', 'w') as f:
    f.write('\n'.join(Seq_nr))

