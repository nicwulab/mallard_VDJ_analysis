import os
import pandas as pd
from Bio.Seq import Seq

def na_remove(tb):
    tb = tb[~tb.fr1.isna()]
    tb = tb[~tb.cdr1.isna()]
    tb = tb[~tb.fr2.isna()]
    tb = tb[~tb.cdr2.isna()]
    tb = tb[~tb.fr3.isna()]
    tb = tb[~tb.cdr3.isna()]
    tb = tb[~tb.fr4.isna()]
    return tb

def FindTheTranslate(input, seq_IG):
    for i in range(3):
        seq = Seq(input[i:]).translate()
        if seq_IG in seq:
            return i, seq

def GetSeqFragment(tb, i):
    input = tb.input.iloc[i]
    fr1 = tb.fr1.iloc[i]
    cdr1 = tb.cdr1.iloc[i]
    fr2 = tb.fr2.iloc[i]
    cdr2 = tb.cdr2.iloc[i]
    fr3 = tb.fr3.iloc[i]
    cdr3 = tb.cdr3.iloc[i]
    fr4 = tb.fr4.iloc[i]
    return input, fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4

def AA2NA(input, seq_aa, fr1):
    ns = seq_aa.find(fr1) * 3 + i_start
    ne = ns + len(fr1) * 3 
    return input[ns:ne]

def PlotFregmentLength(tb, tb_id):
    plt.figure(figsize=(6, 3))
    tb.fr1.str.len().plot(kind='density', label = 'fr1')
    tb.fr2.str.len().plot(kind='density', label = 'fr2') 
    tb.fr3.str.len().plot(kind='density', label = 'fr3')
    tb.fr4.str.len().plot(kind='density', label = 'fr4')
    tb.cdr1.str.len().plot(kind='density', label = 'cdr1')
    tb.cdr2.str.len().plot(kind='density', label = 'cdr2')
    tb.cdr3.str.len().plot(kind='density', label = 'cdr3')
    plt.title(tb_id)
    # set lim for x axis
    plt.xlim(0, 50)
    plt.ylim(0, 5)
    plt.legend()
    plt.savefig(f'plot/fragment_length_{tb_id}.png')

tb_list = os.popen('ls ../2411_4DcukHA/result/*VDJ.tsv').read().split()

for tb_file in tb_list:
    tb = pd.read_csv(tb_file, sep='\t')
    tb = na_remove(tb)
    tb_id = tb_file.split('/')[-1].split('.')[0]
    # plot the length distribution of each fragment
    PlotFregmentLength(tb, tb_id)
    for i in range(len(tb)):
        try:
            seq_id = tb_id + ":" + tb.seq_id.iloc[i]
            input, fr1, cdr1, fr2, cdr2, fr3, cdr3, fr4 = GetSeqFragment(tb, i)
            seq_IG = fr1+cdr1+fr2+cdr2+fr3+cdr3+fr4
            i_start, seq_aa = FindTheTranslate(input, seq_IG)
            with open('result/fr1_3.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, fr1+cdr1+fr2+cdr2+fr3) + '\n')
            with open('result/fr1.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, fr1) + '\n')
            with open('result/cdr1.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, cdr1) + '\n')
            with open('result/fr2.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, fr2) + '\n')
            with open('result/cdr2.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, cdr2) + '\n')
            with open('result/fr3.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, fr3) + '\n')
            with open('result/cdr3.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, cdr3) + '\n')
            with open('result/fr4.fa', 'a') as f:
                f.write(f">{seq_id}\n" + AA2NA(input, seq_aa, fr4) + '\n')
        except:
            pass


fa_list = os.popen('ls result/*.fa').read().split()
for fa in fa_list:
    id = fa.split('/')[-1].split('.')[0]
    # checking if the result is already exist
    if os.path.exists(f"result/{id}.tsv"):
        continue
    cmd = f"blastn -query {fa} -db ../data/20241202_DuckWGS_assemble/Bird75_blast -evalue 1e-4 -max_target_seqs 1 -num_threads 60 -max_hsps 100 -outfmt '6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps' -word_size 15 -qcov_hsp_perc 100 -out result/{id}.tsv"
    os.system(cmd)


