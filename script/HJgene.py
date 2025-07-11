#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
seq_loc="../data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# fr for scaffold 189
TB = pd.read_csv('result/potential_fr_complete.tsv', sep = '\t', index_col = 0)
TB = TB[~TB.fr4_s.isna()]
TB = TB[["189" in i for i in TB.index]]
TB.fr4_s = TB.fr4_s.astype(int)
TB.fr4_e = TB.fr4_e.astype(int)

i = 0
seq_aa = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].translate().seq)
seq_nr = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].seq)
 
seq_nr.find('CGCCGTG') - seq_nr.find('GGTTTTTGG') -len('CGCCGTG')


# signal sequences for D gene
seq_nr_V = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr1_s):TB.iloc[i].fr3_e + 100].seq)
seq_nr_V.find('ACAAAACCC') - seq_nr_V.find('CACGGTG') - len('ACAAAACCC')

Vgene_end = TB.iloc[i].fr3_e + 100 + seq_nr_V.find('ACAAAACCC') + len('ACAAAACCC')
Jgene_start = TB.iloc[i].fr4_s-1-120 + seq_nr.find('CGCCGTG')

import re

Seq_Between = seq_dict[TB.index[i].split('_')[0]][Vgene_end:Jgene_start].seq


indices_1 = [match.start() for match in re.finditer('GATTTTG', str(Seq_Between))]
indices_2 = [match.start() for match in re.finditer('CACCGTG', str(Seq_Between))]
indices_3 = [match.start() for match in re.finditer("CACGGTG", str(Seq_Between))]
indices_4 = [match.start() for match in re.finditer("ACAAAAA", str(Seq_Between))]

[match.start() for match in re.finditer("ACGGT", str(Seq_Between))]


np.array(indices_2) - np.array(indices_1)

indices_1


str(Seq_Between[2083:2161])
str(Seq_Between[3162:3245])
str(Seq_Between[4132:4214]) 
str(Seq_Between[5282:5367]) 
str(Seq_Between[6190:6282])



indices_t = [match.start() for match in re.finditer('TTTT', str(Seq_Between))]
indices_t = indices_1
indices_a = [match.start() for match in re.finditer('AAAAA', str(Seq_Between))]

Paires = []
for i in range(len(indices_t)):
    for ii in range(len(indices_a)):
        if indices_a[ii] - indices_t[i] < 120 and indices_a[ii] - indices_t[i]  > 56:
            # print(indices_t[i]-3, indices_a[ii]+7, indices_a[ii] - indices_t[i] )
            seq_tmp = str(Seq_Between[indices_t[i]-1:indices_a[ii]+7])
            print(indices_t[i], indices_a[ii]+7, seq_tmp[:9], seq_tmp[9+12:9+12+7], 
                  len(seq_tmp) - 2*(9+12+7), seq_tmp[-9-12-6:-9-12], seq_tmp[-9:])
            Paires.append([indices_t[i]-3, indices_a[ii]+7, seq_tmp[:9], seq_tmp[9+12:9+12+7], 
                  len(seq_tmp) - 2*(9+12+7), seq_tmp[-9-12-6:-9-12], seq_tmp[-9:]])


Seq_Vgene = seq_dict[TB.index[0].split('_')[0]][:Vgene_end].seq

for pair in Paires:
    sig1 = Seq(pair[3]).reverse_complement()
    sig2 = Seq(pair[2]).reverse_complement()
    print(sig1, sig2, pair[3], pair[2])

    # sig3 = pair[5]
    # sig4 = pair[6]
    indices_1 = [match.start() for match in re.finditer(str(sig1), str(Seq_Vgene))]
    indices_2 = [match.start() for match in re.finditer(str(sig2), str(Seq_Vgene))]
    if len(indices_1) != 0 and len(indices_2) != 0:
        # pariwise distance between the two signals
            for i in indices_1:
                for ii in indices_2:
                    if ii -i >0 and ii - i < 500:
                        print(i, ii, ii - i)










## Check the D gene in Chicken
seq_loc="data/IMGT000014_FASTA_20250314_22h27m34_2458062731166920083.txt"
seq2_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))
Chicken = list(seq2_dict.items())[0][1].seq

Chicken.find('gtactgctggtagcatcgacgcatggggccacgggaccgaagtcatcgtctcctccg')

str(Chicken[95048:95348].translate())
Seq_Between = Chicken[95348:106171]

# 95348-95391, 106171-106133

indices_1 = [match.start() for match in re.finditer('gattttg', str(Seq_Between))]
indices_2 = [match.start() for match in re.finditer('caccgtg', str(Seq_Between))]
indices_3 = [match.start() for match in re.finditer("cacggtg", str(Seq_Between))]
indices_4 = [match.start() for match in re.finditer("acaaaaa", str(Seq_Between))]

[i + 95348 for i in indices_1]

indices_t = [match.start() for match in re.finditer('tttt', str(Seq_Between))]
indices_t = indices_1
indices_a = [match.start() for match in re.finditer('aaaaa', str(Seq_Between))]

for i in range(len(indices_t)):
    for ii in range(len(indices_a)):
        if indices_a[ii] - indices_t[i] < 120 and indices_a[ii] - indices_t[i]  > 56:
            # print(indices_t[i]-3, indices_a[ii]+7, indices_a[ii] - indices_t[i] )
            seq_tmp = str(Seq_Between[indices_t[i]0:indices_a[ii]+7])
            print(indices_t[i] + 95348, indices_a[ii]+7 + 95348, seq_tmp[:9],
                  seq_tmp[9+12:9+12+7], len(seq_tmp) - 2*(9+12+7), seq_tmp[-9-12-6:-9-12], seq_tmp[-9:])


len(Seq_Between)


tmp = Seq('ggattttgg')
tmp = tmp.reverse_complement()
tmp2 = Seq('caccgtg')
tmp2 = tmp2.reverse_complement()
indices_all1 = [match.start() for match in re.finditer(str(tmp), str(Chicken))]
indices_all2 = [match.start() for match in re.finditer(str(tmp2), str(Chicken))]

tmp = Seq('cgattttgg')
tmp = tmp.reverse_complement()
tmp2 = Seq('gattcta')
tmp2 = tmp2.reverse_complement()
indices_all1 = [match.start() for match in re.finditer(str(tmp), str(Chicken))]
indices_all2 = [match.start() for match in re.finditer(str(tmp2), str(Chicken))]


# fr for scaffold 063 Light chain
TB = pd.read_csv('result/potential_fr_complete.tsv', sep = '\t', index_col = 0)
TB = TB[~TB.fr4_s.isna()]
TB = TB[["063" in i for i in TB.index]]
TB.fr4_s = TB.fr4_s.astype(int)
TB.fr4_e = TB.fr4_e.astype(int)

i = 0

str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].translate().seq)

seq_aa = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].translate().seq)
seq_nr = str(seq_dict[TB.index[i].split('_')[0]][(TB.iloc[i].fr4_s-1-120):TB.iloc[i].fr4_e].seq)
 
seq_nr.find('CGCCGTG') - seq_nr.find('GGTTTTTGG')










# Check the J gene by the way
import re

Potential_all = pd.read_csv('result/potential_fr_complete.tsv', sep = '\t', index_col = 0)
FW4 = Potential_all[~Potential_all.fr4_s.isna()]
FW4 = FW4[FW4.index == 'ptg000189l_17'].iloc[0]
FW4.fr4_s = FW4.fr4_s.astype(int)
FW4.fr4_e = FW4.fr4_e.astype(int)

J_ge = str(seq_dict[FW4.name.split('_')[0]][FW4.fr4_s-301:FW4.fr4_e].seq)
J_rc = str(seq_dict[FW4.name.split('_')[0]][FW4.fr4_s-301:FW4.fr4_e].seq.reverse_complement())

js = J_ge.find('GGTTTTTGG')
Je = J_ge.find('CGCCGTG')

JS = seq_dict[FW4.name.split('_')[0]].seq.find(J_ge[js:])
JE = JS + len(J_ge[js:]) 
J_gene = str(seq_dict[FW4.name.split('_')[0]].seq[JS:JE])
