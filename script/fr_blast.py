#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from collections import Counter

fr_list = os.popen('ls result/fr*.tsv').read().split()

fr_all = pd.DataFrame()
for fr_file in fr_list:
    fr = fr_file.split('/')[-1].split('.')[0]
    tb = pd.read_csv(fr_file, sep='\t')
    tb.columns = ['qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps'] 
    tb['seg'] = fr
    fr_all = pd.concat([fr_all, tb])

# split by scalfold and identify the directions
DictV = {}
scalfolds = fr_all[fr_all.seg == 'fr4'].sacc.unique()
for sfd in scalfolds:
    tb = fr_all[fr_all.sacc == sfd]
    tb['forward'] = tb.sstart[tb.seg == 'fr1'].mean() - tb.sstart[tb.seg == 'fr4'].mean() < 0
    tb['seg_f'] = tb.sstart - tb.send < 0
    # remove the reversed 
    tb = tb[tb.forward == tb.seg_f]
    # remove duplicates 
    if tb.forward.iloc[0]:
        tb = tb.sort_values(['sstart', 'send'], ascending=[True, False])
    else:
        tb = tb.sort_values(['sstart', 'send'], ascending=[False, True])
    # remove duplicates
    tb = tb[~tb.sstart.duplicated()]
    tb = tb[~tb.send.duplicated()]
    # sort again 
    if tb.forward.iloc[0]:
        tb = tb.sort_values(['sstart', 'send'], ascending=[True, False])
    else:
        tb = tb.sort_values(['sstart', 'send'], ascending=[False, True])
    # remove overlaps
    if tb.forward.iloc[0]:
        mask_overlap = np.where(tb.sstart[1:].to_numpy() - tb.send[:-1].to_numpy() < 0)[0]
    else:
        mask_overlap = np.where(tb.sstart[1:].to_numpy() - tb.send[:-1].to_numpy() > 0)[0]
    mask = [i for i in range(tb.shape[0]) if i not in 1+mask_overlap]
    tb = tb.iloc[mask, :]
    # itteration by the frame work
    N = 0
    sacc = tb.sacc.iloc[0]
    DictV[f"{sacc}_{0}"] = {'forward': tb.forward.iloc[0]}
    for i in range(tb.shape[0]):
        seg = tb.seg.iloc[i]
        print(seg)
        if seg == 'fr1':
            N += 1
            DictV[f"{sacc}_{N}"] = {'forward': tb.forward.iloc[i]}
        DictV[f"{sacc}_{N}"][f"{seg}_s"] = tb.sstart.iloc[i]
        DictV[f"{sacc}_{N}"][f"{seg}_e"] = tb.send.iloc[i]


TB_all = pd.DataFrame(DictV).T
TB_all = TB_all[TB_all.isna().sum(axis=1)<8]
TB_all.to_csv('result/potential_fr_all.tsv', sep='\t')
TB_complete = TB_all[TB_all[['fr1_s', 'fr1_e', 'fr2_s', 'fr2_e', 'fr3_s', 'fr3_e']].isna().sum(axis = 1) == 0]
TB_complete['fr_dis'] = TB_complete.fr3_e - TB_complete.fr1_s
TB_complete.fr_dis[TB_complete.forward == False] *= -1    
TB_complete = TB_complete[TB_complete.fr_dis<=320]
TB_complete['aa_fw123'] = '' 
TB_complete['nr_fw123'] = '' 
TB_complete['nr_tail'] = '' 

## take the seq from the genome
from Bio import SeqIO
seq_loc="../data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

for i in range(73):
    if TB_complete.forward.iloc[i]:
        TB_complete.aa_fw123.iloc[i] = str(seq_dict[TB_complete.index[i].split('_')[0]][(TB_complete.iloc[i].fr1_s-1):TB_complete.iloc[i].fr3_e].translate().seq)
        TB_complete.nr_fw123.iloc[i] = str(seq_dict[TB_complete.index[i].split('_')[0]][(TB_complete.iloc[i].fr1_s-1):TB_complete.iloc[i].fr3_e].seq)
        TB_complete.nr_tail.iloc[i] = str(seq_dict[TB_complete.index[i].split('_')[0]][TB_complete.iloc[i].fr3_e:TB_complete.iloc[i].fr3_e+100].seq)
    else:
        TB_complete.aa_fw123.iloc[i] = str(seq_dict[TB_complete.index[i].split('_')[0]][TB_complete.iloc[i].fr3_e-1: TB_complete.iloc[i].fr1_s].reverse_complement().translate().seq)
        TB_complete.nr_fw123.iloc[i] = str(seq_dict[TB_complete.index[i].split('_')[0]][TB_complete.iloc[i].fr3_e-1: TB_complete.iloc[i].fr1_s].reverse_complement().seq)
        TB_complete.nr_tail.iloc[i] = str(seq_dict[TB_complete.index[i].split('_')[0]][TB_complete.iloc[i].fr1_s : TB_complete.iloc[i].fr1_s+100].reverse_complement().seq)

TB_complete.to_csv('result/potential_fr_complete.tsv', sep='\t') 
