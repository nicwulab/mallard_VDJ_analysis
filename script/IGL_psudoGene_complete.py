#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 17 2025

This script is used to extract the IGH pseudogene sequences from the IGH gene
The previous methoud which used to defeine the pseudogene is not accurate.
Partial head of them are missing.
Tail of the pseudogene are not defined.

For the head completenese, we could use the Functional gene to define.
For the tail, we could do some test on D gene and V gene. 
"""

import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio import AlignIO


def p_align(seq1, seq2, gap_open = -10, gap_extend = -0.5):
    matrix = substitution_matrices.load("BLOSUM62")
    alignments = pairwise2.align.localds(seq1, seq2, matrix, gap_open, gap_extend)
    return alignments[0]

def StringCut(my_string):
    # Step 1: split into lines of 60 chars
    # Step 2: split each line into groups of 3 chars
    lines = [my_string[i:i+60] for i in range(0, len(my_string), 60)]
    formatted_lines = [' '.join([line[j:j+3] for j in range(0, len(line), 3)]) for line in lines]
    # Combine back into a single string
    # final_result = '\n'.join(formatted_lines)
    return formatted_lines

def AlignSort(alignAA, NT):
    seq_raw = NT.replace('-', '')
    seq_aan1 = str(Seq(seq_raw[2:]).translate()) + '--'
    seq_aap1 = str(Seq(seq_raw[1:]).translate()) + '--'
    SAP1 = ""
    SAN1 = "-"
    N = 0
    for i in alignAA:
        if i == '-':
            SAP1 += '-'
            SAN1 += '-'
        else:
            SAP1 += seq_aap1[N]
            SAN1 += seq_aan1[N]
            N += 1
    N = 0
    seqA_nt = '' 
    seqA_aa = ''
    seqA_aaN1 = ''
    seqA_aaP1 = ''
    ii = -1
    for i in alignAA:
        ii += 1
        if i == '-':
            seqA_nt += '---'
            seqA_aa += f'-{i}-'
            seqA_aaN1 += f'---' 
            seqA_aaP1 += f'---' 
        else:
            seqA_nt += NT[N:N+3]
            seqA_aa += f' {i} ' 
            seqA_aaN1 += f'{SAN1[ii]}  '
            seqA_aaP1 += f'  {SAP1[ii]}'
            N += 3
    return StringCut(seqA_nt), StringCut(seqA_aa), StringCut(seqA_aaN1), StringCut(seqA_aaP1)

def RSS_Check():
    # RSS Check
    # blastn -query data/RSS_IGL.fa  -db /data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/ptg000063l  -evalue 30000 -outfmt '6 qacc sacc pident qcovs qlen qstart qend sstart send qseq sseq' -word_size 4 > result/RSS_IGL.tsv
    RSS_IGL = pd.read_csv('result/RSS_IGL.tsv', sep='\t', header=None)
    RSS_IGL.columns = ['qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq']
    # remove gap/insertion 
    RSS_IGL = RSS_IGL[~RSS_IGL.sseq.str.contains('-')]
    RSS_IGL = RSS_IGL[~RSS_IGL.qseq.str.contains('-')]
    # add direction
    RSS_IGL['dir'] = np.where(RSS_IGL.sstart < RSS_IGL.send, '+', '-')

    # adjust the star and end base on the qstart and qend
    RSS_IGL.sstart[RSS_IGL.dir =='+'] -=RSS_IGL.qstart[RSS_IGL.dir =='+'] - 1
    RSS_IGL.send[RSS_IGL.dir =='+'] += RSS_IGL.qlen[RSS_IGL.dir =='+'] - RSS_IGL.qend[RSS_IGL.dir =='+']

    RSS_IGL.sstart[RSS_IGL.dir =='-'] +=RSS_IGL.qstart[RSS_IGL.dir =='-'] - 1
    RSS_IGL.send[RSS_IGL.dir =='-'] -= RSS_IGL.qlen[RSS_IGL.dir =='-'] - RSS_IGL.qend[RSS_IGL.dir =='-']

    # (RSS_IGL.sstart - RSS_IGL.send).unique()
    RSS_IGL1 = RSS_IGL[RSS_IGL.qacc =='RSS_IGL1']
    RSS_IGL2 = RSS_IGL[RSS_IGL.qacc =='RSS_IGL2']

    RSS_IGL1_fend = RSS_IGL1.send[RSS_IGL1.dir == '+']
    RSS_IGL1_rstart = RSS_IGL1.sstart[RSS_IGL1.dir == '-']
    RSS_IGL2_fstart = RSS_IGL2.sstart[RSS_IGL2.dir == '+']
    RSS_IGL2_rend = RSS_IGL2.send[RSS_IGL2.dir == '-']

    # Pairwise distance check
    RSSPair_f = pd.DataFrame([RSS_IGL1_fend.to_list() + RSS_IGL2_fstart.to_list(),["RSS1"]*len(RSS_IGL1_fend)+ ["RSS2"] * len(RSS_IGL2_fstart)], index = ['pos', 'Type']).T
    RSSPair_f.sort_values('pos', inplace=True, ignore_index=True)
    RSSPair_f.pos = RSSPair_f.pos.astype(int)

    Mask = [i for i in range(RSSPair_f.shape[0]-1) if RSSPair_f.Type[i] == 'RSS1' and RSSPair_f.Type[i+1] == 'RSS2']
    RSSPair_f_adj = RSSPair_f.iloc[Mask, :]
    RSSPair_f_adj2 = RSSPair_f.iloc[[i+1 for i in Mask], :]
    RSS_diff_f = RSSPair_f_adj2.pos.to_numpy() - RSSPair_f_adj.pos.to_numpy()

    tmp22 = RSSPair_f_adj.iloc[np.where(RSS_diff_f==23)[0],:]
    tmp23 = RSSPair_f_adj.iloc[np.where(RSS_diff_f==24)[0],:]
    tmp24 = RSSPair_f_adj.iloc[np.where(RSS_diff_f==25)[0],:]
    tmp22['Spacer'] = 22
    tmp23['Spacer'] = 23
    tmp24['Spacer'] = 24
    RSSPair_f_result = pd.concat([tmp22, tmp23, tmp24], axis=0)
    RSSPair_f_result.sort_values('pos', inplace=True, ignore_index=True)
    RSSPair_f_result.pos -= 7
    return RSSPair_f_result


#####################
## Preparing the data 
#####################

# read the genome
seq_loc="/data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# preparing the IGL gene
CMD = '''awk '{print ">"FILENAME$2"\\n"$3}' ../2411_4DcukHA/result/*VDJ.tsv |sed 's=../2411_4DcukHA/result/==;s/.tsv/:/' > result/fullVDJ.fasta'''
os.system(CMD)

seq_loc="result/fullVDJ.fasta"
fVDJSeq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# blast into to chrom
CMD = '''makeblastdb -parse_seqids -dbtype nucl -in data/2chrom.fa -out blastdb/2chrom'''
os.system(CMD)

# blast all seq to 2 scaffolds  
CMD = '''blastn -query result/fullVDJ.fasta  -db blastdb/2chrom -evalue 1e-5 -max_target_seqs 1 -num_threads 60 -max_hsps 1 -outfmt '6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps' > result/blast2chrom.tsv'''
os.system(CMD)

# split to IGH and IGV
tb_hl = pd.read_csv('result/blast2chrom.tsv', sep='\t')
tb_hl.columns = ['qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps']

with open('result/IGH.list', 'w') as f:
    f.write("\n".join(tb_hl.qacc[tb_hl.sacc == 'ptg000189l'].to_list()))
with open('result/IGL.list', 'w') as f:
    f.write("\n".join(tb_hl.qacc[tb_hl.sacc == 'ptg000063l'].to_list()))

os.system('seqkit grep -f result/IGH.list result/fullVDJ.fasta > result/IGH.fa')
os.system('seqkit grep -f result/IGL.list result/fullVDJ.fasta > result/IGL.fa')

seq_loc="result/IGL.fa"
VDJL_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))


#####################
## J gene
#####################

chrm = "ptg000063l"

RSS_IGLJ = 'GGTTTTTGCATGGCCCTGTATCACTGTGAGGTGTATTTGGGGCCGGGACCACGTTGACCGTCCTG'
IGLJ = 'AGGTGTATTTGGGGCCGGGACCACGTTGACCGTCCTG'

J_start = seq_dict[chrm].seq.find(IGLJ)
J_end = J_start + len(IGLJ)
with open('final/IGLJ.fa', 'w') as f:
    f.write(f">IGLJ1 {chrm} {J_start+1} {J_end}\n{IGLJ}\n")

seq_loc="final/IGLJ.fa"
IGLJ_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# D gene annotation
## make the blastdb 
CMD = 'makeblastdb -parse_seqids -dbtype nucl -in final/IGLJ.fa  -out blastdb/IGLJ'
os.system(CMD)

## balst the VDJ-Seq to the D genes
CMD = "blastn -strand plus -query result/IGL.fa  -db blastdb/IGLJ -evalue 1 -max_target_seqs 1 -num_threads 60 -max_hsps 1  -word_size 4 -outfmt '6 qacc sacc pident evalue qcovs qlen qstart qend sstart send length mismatch gaps' -out result/IGL-IGLJ.tsv"
os.system(CMD)

tb_vdj_j = pd.read_csv('result/IGL-IGLJ.tsv', sep='\t', header=None)
tb_vdj_j.columns = ['qacc', 'sacc', 'pident', 'evalue', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps']

Leader_IGLV = []
for i in range(tb_vdj_j.shape[0]):
    seq_tmp = tb_vdj_j.iloc[i]
    seq = str(VDJL_dict[seq_tmp.qacc].seq)[:seq_tmp.qstart - seq_tmp.sstart]
    Leader_IGLV.append(f">{seq_tmp.qacc}\n{seq}")

with open('result/Leader_IGLV.fa', 'w') as f:
    f.write('\n'.join(Leader_IGLV))
#####################
## V gene
#####################

# Functional V genes
seq_loc = 'result/IGLV_fun.fa'
IGLV_fun_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

[seq_dict['ptg000063l'].seq.find(seq.seq) for seq in IGLV_fun_dict.values()]

IGLV1 = 'CAGGCAGCGCTGACTCAGCCGGCCTCGAAGTCGGTGAATCCGGGAGACACCGTGCAGATCACTTGCTCCGGGGGTGGCAGCTACTACGGCTGGTTCCAGCAGAAGACCCCTGGCACTGGCCCTGTCACCGTGATCTATGACAATACCAACAGACCCTCGGGCATCCCTTCTCGATTCTCTGCTTCCACATCTGGCTCCGTGTCCACTTTAACCATCACTGGGGTCCAAGCCGAGGACGAGGCTGTCTATTACTGTGGTGACTACAGTAG'

# Cut the Leader based on the Functional V gene
# 1. blast the Functional V gene to the potential V gene
## Make the blastdb
CMD = 'makeblastdb -parse_seqids -dbtype nucl -in result/IGLV_fun.fa  -out blastdb/IGLV_fun'
os.system(CMD)
## blast
CMD = "blastn -strand plus -query result/Leader_IGLV.fa  -db blastdb/IGLV_fun -evalue 1 -max_target_seqs 100000 -num_threads 60 -max_hsps 10000  -word_size 4 -outfmt '6 qacc sacc pident evalue qcovs qlen qstart qend sstart send length mismatch gaps' -out result/LIGLV-IGLV_fun.tsv -qcov_hsp_perc 70"
os.system(CMD)

# 2. Find the head of the V gene
tb_vdj_v = pd.read_csv('result/LIGLV-IGLV_fun.tsv', sep='\t', header=None)
tb_vdj_v.columns = ['qacc', 'sacc', 'pident', 'evalue', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps']

# 3. save the V gene
seq_loc="result/Leader_IGLV.fa"
Leader_H_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))
IGLV = []
for i in range(tb_vdj_v.shape[0]):
    qacc = tb_vdj_v.qacc.iloc[i]
    sacc = tb_vdj_v.sacc.iloc[i]
    sstart = tb_vdj_v.sstart.iloc[i]
    qstart = tb_vdj_v.qstart.iloc[i]
    if sstart >= qstart:
        seq = Leader_H_dict[qacc].seq[qstart-sstart:]
        if len(seq) >= 200:
            IGLV.append(f">{sacc}\n{seq}")


with open('result/IGLV.fa', 'w') as f:
    f.write('\n'.join(IGLV))

# 4. align the V gene back to the scaffolds
CMD = 'blastn -query result/IGLV.fa  -db /data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/ptg000063l -evalue 1e-5 -max_target_seqs 100 -num_threads 1 -max_hsps 100 -outfmt "6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps sseq" -qcov_hsp_perc 70 -out result/IGLV-ptg000063l.tsv'
os.system(CMD)

# 5. read the blast result
tb_vdj_genome = pd.read_csv('result/IGLV-ptg000063l.tsv', sep='\t', header=None)
tb_vdj_genome.columns = ['qacc', 'sacc', 'pident', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps', 'sseq']
tb_vdj_genome['rela_pos'] =  ((tb_vdj_genome.sstart +tb_vdj_genome.send )/2/10).astype(int)
tb_vdj_genome['MAX'] = tb_vdj_genome[['sstart', 'send']].max(axis=1)
tb_vdj_genome['MIN'] = tb_vdj_genome[['sstart', 'send']].min(axis=1)
tb_vdj_genome.sort_values('MIN', ascending=False, inplace=True)

Switch = np.concatenate([[False], tb_vdj_genome.MAX.to_numpy()[1:] < tb_vdj_genome.MIN.to_numpy()[:-1]])
Slist = np.where(Switch==True)[0]
Slist = np.concatenate([[0], Slist])

tb_vdj_genome.index = range(tb_vdj_genome.shape[0])

geneType = {0:'+ 1',
            1:'+ 2',
            2:'+ 3',
            3:'- 1',
            4:'- 2',
            5:'- 3',}

def Align_find(align_best, seq_nt):
    tmp1, tmp2, _, _ = AlignSort(align_best.seqA, IGLV1)
    tmp3, tmp4,tmp5,tmp6 = AlignSort(align_best.seqB, seq_nt)
    Result = "\n".join(['\n'.join([tmp6[i],tmp5[i],tmp4[i],tmp2[i], tmp1[i], tmp3[i], '\n']) for i in range(len(tmp3))])
    return Result

from collections import Counter

def Align_best(MIN, MAX):
    seq1 = str(seq_dict[chrom].seq[MIN:MAX].translate())
    seq2 = str(seq_dict[chrom].seq[MIN+1:MAX].translate())
    seq3 = str(seq_dict[chrom].seq[MIN+2:MAX].translate())
    seq4 = str(seq_dict[chrom].seq[MIN:MAX].reverse_complement().translate())
    seq5 = str(seq_dict[chrom].seq[MIN:MAX-1].reverse_complement().translate())
    seq6 = str(seq_dict[chrom].seq[MIN:MAX-2].reverse_complement().translate())
    seq0 = str(Seq(IGLV1).translate())
    SeqAll = [seq1, seq2, seq3, seq4, seq5, seq6]
    align_all = [p_align(seq0, seq2, gap_open = -10, gap_extend = -.5) for seq2 in SeqAll]
    score_swith = np.array([algin.score for algin in align_all])
    mask_max = np.where(score_swith == max(score_swith))[0][0]
    direct = geneType[mask_max].split()[0]
    # start_tmp = 0 
    jump = int(geneType[mask_max].split()[-1])
    align_best = align_all[mask_max]
    head_gap = 0
    for i in align_best.seqB:
        if i != '-':
            break
        head_gap += 1
    return direct, jump, head_gap, align_best


RSSPair_f_result = RSS_Check()
chrom = tb_vdj_genome.sacc.iloc[0]
Seq_aa = []
Seq_nt = []
N = 0
for i in range(Slist.shape[0]-1):
    tmp = tb_vdj_genome.iloc[Slist[i]:Slist[i+1], :]
    # check fo 100 qcovs first and then
    tmp = tmp[tmp.qcovs>=90]
    # take the longest one
    MAX = tmp.MAX.max()
    MIN = tmp.MIN.min()
    # best alignment result
    direct, jump, head_gap, align_best = Align_best(MIN, MAX)
    if align_best.score> 50:
        # fill the gap in the head
        if direct == '+':
            MIN +=  jump -1
            direct, jump, head_gap, align_best = Align_best(MIN, MAX)
            MIN -= head_gap  * 3 
            direct, jump, head_gap, align_best = Align_best(MIN, MAX)
        if direct == '-':
            MAX -=  jump -1
            direct, jump, head_gap, align_best = Align_best(MIN, MAX)
            MAX += head_gap  * 3 
            direct, jump, head_gap, align_best = Align_best(MIN, MAX)
        N += 1
        id = f">IGLV1-{N} {chrom} {MIN} {MAX} {direct}"
        # start_tmp = 0 
        if '+' == direct:
            seq_nt = str(seq_dict[chrom].seq[MIN :MAX])
            seq_aa = align_best.seqB.replace('-', '')
            # Rss Check
            RSS_cloest = min(abs(RSSPair_f_result.pos - MAX))
            if RSS_cloest < 30:
                RSS_tmp = RSSPair_f_result[abs(RSSPair_f_result.pos - MAX)==RSS_cloest].iloc[0, :].to_dict()
                RSS_seq = seq_dict[chrom].seq[RSS_tmp['pos']:RSS_tmp['pos']+RSS_tmp['Spacer']+7+9]
                RSS_Spacer = RSS_tmp['Spacer']
                id = f"{id} {RSS_seq} {RSS_Spacer} {RSS_cloest}"
                direct, jump, head_gap, align_best = Align_best(MIN, RSS_tmp['pos'])
                seq_nt = str(seq_dict[chrom].seq[MIN :RSS_tmp['pos']])
                seq_aa = align_best.seqB.replace('-', '')
        else:
            seq_nt = str(seq_dict[chrom].seq[MIN:MAX].reverse_complement())
            seq_aa = align_best.seqB.replace('-', '')
        Result = Align_find(align_best, seq_nt)
        Seq_aa.append(f"{id}\n{seq_aa}")
        Seq_nt.append(f"{id}\n{seq_nt}")
        with open(f'result/IGLV_align/IGLV1-{N}.aln', 'w') as f:
            f.write(id+'\n')
            f.write(Result)

with open('final/IGLV.fa', 'w') as f:
    f.write('\n'.join(Seq_nt))
with open('final/IGLV.pro', 'w') as f:
    f.write('\n'.join(Seq_aa))

# muscle align the V gene

CMD = 'muscle -align final/IGLV.fa -output result/IGLV.aln'
os.system(CMD)

CMD = 'muscle -align final/IGLV.pro -output result/IGLV.pro.aln'
os.system(CMD)
