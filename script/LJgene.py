#!/usr/bin/env python3
import pandas as pd
import re

from Bio import SeqIO
seq_loc="../data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# fr for scaffold 189
TB = pd.read_csv('result/potential_fr_complete.tsv', sep = '\t', index_col = 0)
TB = TB[~TB.fr4_s.isna()]
TB = TB[["63l" in i for i in TB.index]]
TB.fr4_s = TB.fr4_s.astype(int)
TB.fr4_e = TB.fr4_e.astype(int)

i = 0
seq_aa = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].translate().seq)
seq_nr = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].seq)
 
NT = seq_nr.find('TTTTT')

seq_nr[NT-2:NT+7+23+10]

len(seq_nr[NT-2:])

seq_dict[TB.index[i].split('_')[0]].seq.find(seq_nr[NT-2:])

seq_dict[TB.index[i].split('_')[0]][9508724:9508789].seq.translate()


# check them on the V gene
import os
import numpy as np
import pandas as pd
from collections import Counter

fr_list = [i for i in os.popen('ls result/ptg00*.tsv').read().split() if 'fr1_3'  in i][:1]

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


for i in range(30): #tb.shape[0]):
    tmp = tb.iloc[i]
    seq = seq_dict[tmp.sacc]
    if tmp.sstart < tmp.send:
        # print("FW", seq[tmp.sstart-1:tmp.send].translate().seq)
        print("FW", seq[tmp.sstart-1:tmp.send].seq)
    else:
        # print("RW", seq[tmp.send-1:tmp.sstart+51+12].reverse_complement().translate().seq)
        print("RW", seq[tmp.send-1:tmp.sstart].reverse_complement().seq)
        # print("RW", seq[tmp.send-1:tmp.sstart+51+12].reverse_complement().seq)


    seq[tmp.sstart:tmp.send].seq.translate()



