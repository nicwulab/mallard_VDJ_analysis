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

D_gene= ['AGTGCTGGTGCTGGTTGGGTTA',
         'GGTACTGGTGGTTGGGGTGCTAATGCA',
         'AGTACTGGTGCTTGTTGTGGTGGTTA',
         'GGTGCTGGTAGTGGTTACGGTGCTTATTA',
         'GGTTATGCTAGTTGTGGTGGTTATACTTGTGCTTAT']

IGHV1 = 'GCTGCCACCTTGGATGAGTCCGGAGGGGGCCTCGTGAGTCCCGGGGGGTCCCTGACCCTGGTCTGCAAGGGCTCCGGATACACCTTCAGCAGCTACGGCATGGGATGGGTGCGACAGGCACCCGGGAAGGGGCTCGAGTACGTCGCGAGTATTAACAGCAGTGGTAGTAGCACTTACTACGCGCCGGCGGTGAAGGGACGCTTCACCATCTCCAGGAACAACGGGCAGAGCACGCTCACCCTGCAGATGAACAGCCTCAAGGCCGAAGACACCGCCACCTACTACTGCGCGAAAGCTGCTGGT'

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



#####################
## Preparing the data 
#####################

# read the genome
seq_loc="../data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"
seq_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# preparing the IGH gene
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

seq_loc="result/IGH.fa"
VDJH_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))


#####################
## D gene
#####################

# [(2083,2160), AGTGCTGGTGCTGGTTGGGTTA),
#  (3162,3244), GGTACTGGTGGTTGGGGTGCTAATGCA), 
#  (4132,4213), AGTACTGGTGCTTGTTGTGGTGGTTA), 
#  (5282,5366), GGTGCTGGTAGTGGTTACGGTGCTTATTA), 
#  (6190,6281), GGTTATGCTAGTTGTGGTGGTTATACTTGTGCTTAT)]

List_D = [(2083,2160), (3162,3244),  (4132,4213), (5282,5366), (6190,6281)]

with open('final/IGHD.fa', 'w') as f:
    N = 0
    for i in List_D:
        N += 1
        dgene = str(seq_dict[fun_gene.sacc].seq[i[0]+534855 +12+9+7-1: i[1]+534855-9-12-7])
        f.write(f">IGHD{N} {tb_blast.iloc[0].sacc} {i[0]} {i[1]} \n{dgene}\n")

seq_loc="final/IGHD.fa"
IGHD_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))

# D gene annotation
## make the blastdb 
CMD = 'makeblastdb -parse_seqids -dbtype nucl -in final/IGHD.fa  -out blastdb/IGHD'
os.system(CMD)

## balst the VDJ-Seq to the D genes
CMD = "blastn -strand plus -query result/IGH.fa  -db blastdb/IGHD -evalue 1 -max_target_seqs 1 -num_threads 60 -max_hsps 1  -word_size 4 -outfmt '6 qacc sacc pident evalue qcovs qlen qstart qend sstart send length mismatch gaps' -out result/IGH-IGHD.tsv"
os.system(CMD)

tb_vdj_d = pd.read_csv('result/IGH-IGHD.tsv', sep='\t', header=None)
tb_vdj_d.columns = ['qacc', 'sacc', 'pident', 'evalue', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps']

Leader_IGHV = []
for i in range(tb_vdj_d.shape[0]):
    seq_tmp = tb_vdj_d.iloc[i]
    seq = str(VDJH_dict[seq_tmp.qacc].seq)[:seq_tmp.qstart - seq_tmp.sstart]
    Leader_IGHV.append(f">{seq_tmp.qacc}\n{seq}")

with open('result/Leader_IGHV.fa', 'w') as f:
    f.write('\n'.join(Leader_IGHV))


#####################
## V gene
#####################

with open('result/IGHV1-1_99.fa', 'w') as f:
    f.write(f">IGHV1-1\n{IGHV1[:99]}\n")

# Cut the Leader based on the Functional V gene
# 1. blast the Functional V gene to the potential V gene
## Make the blastdb
CMD = 'makeblastdb -parse_seqids -dbtype nucl -in result/Leader_IGHV.fa  -out blastdb/Leader_IGHV'
os.system(CMD)
## blast
CMD = "blastn -strand plus -query result/IGHV1-1_99.fa  -db blastdb/Leader_IGHV -evalue 1 -max_target_seqs 100000 -num_threads 60 -max_hsps 10000  -word_size 4 -outfmt '6 qacc sacc pident evalue qcovs qlen qstart qend sstart send length mismatch gaps' -out result/IGHV1-1_99-Leader_IGHV.tsv -qcov_hsp_perc 70"
os.system(CMD)

# 2. Find the head of the V gene
tb_vdj_v = pd.read_csv('result/IGHV1-1_99-Leader_IGHV.tsv', sep='\t', header=None)
tb_vdj_v.columns = ['qacc', 'sacc', 'pident', 'evalue', 'qcovs', 'qlen', 'qstart', 'qend', 'sstart', 'send', 'length', 'mismatch', 'gaps']

# 3. save the V gene
seq_loc="result/Leader_IGHV.fa"
Leader_H_dict = SeqIO.to_dict(SeqIO.parse(seq_loc, "fasta"))
IGHV = []
for i in range(tb_vdj_v.shape[0]):
    sacc = tb_vdj_v.sacc.iloc[i]
    sstart = tb_vdj_v.sstart.iloc[i]
    qstart = tb_vdj_v.qstart.iloc[i]
    if sstart >= qstart:
        seq = Leader_H_dict[sacc].seq[sstart-qstart:]
        if len(seq) >= 200:
            IGHV.append(f">{sacc}\n{seq}")
with open('result/IGHV.fa', 'w') as f:
    f.write('\n'.join(IGHV))

# 4. align the V gene back to the scaffolds
CMD = 'blastn -query result/IGHV.fa  -db ../data/20241202_DuckWGS_assemble/ptg000189l -evalue 1e-5 -max_target_seqs 100 -num_threads 1 -max_hsps 100 -outfmt "6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps sseq" -qcov_hsp_perc 70 -out result/IGHV-ptg000189l.tsv'
os.system(CMD)

# 5. read the blast result
tb_vdj_genome = pd.read_csv('result/IGHV-ptg000189l.tsv', sep='\t', header=None)
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
    tmp1, tmp2, _, _ = AlignSort(align_best.seqA, IGHV1)
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
    seq0 = str(Seq(IGHV1).translate())
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


chrom = tb_vdj_genome.sacc.iloc[0]
Seq_aa = []
Seq_nt = []
N = 0
for i in range(58): #Slist.shape[0]-1):
    tmp = tb_vdj_genome.iloc[Slist[i]:Slist[i+1], :]
    # check fo 100 qcovs first and then
    tmp = tmp[tmp.qcovs>=100]
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
        id = f">IGHV1-{N} {chrom} {MIN} {MAX} {direct}"
        # start_tmp = 0 
        if '+' == direct:
            seq_nt = str(seq_dict[chrom].seq[MIN :MAX])
            seq_aa = align_best.seqB.replace('-', '')
        else:
            seq_nt = str(seq_dict[chrom].seq[MIN:MAX].reverse_complement())
            seq_aa = align_best.seqB.replace('-', '')
        Result = Align_find(align_best, seq_nt)
        Seq_aa.append(f"{id}\n{seq_aa}")
        Seq_nt.append(f"{id}\n{seq_nt}")
        with open(f'result/IGHV_align/IGHV1-{N}.aln', 'w') as f:
            f.write(id+'\n')
            f.write(Result)

with open('final/IGHV.fa', 'w') as f:
    f.write('\n'.join(Seq_nt))
with open('final/IGHV.pro', 'w') as f:
    f.write('\n'.join(Seq_aa))

# muscle align the V gene

CMD = 'muscle -align final/IGHV.fa -output result/IGHV.aln'
os.system(CMD)

CMD = 'muscle -align final/IGHV.pro -output result/IGHV.pro.aln'
os.system(CMD)


#####################
## J gene
#####################

sstart = 541834
send = 541928

jgene = str(seq_dict[fun_gene.sacc].seq[sstart+7+9+23-1: send])

with open('final/IGHJ.fa', 'w') as f:
    f.write(f">IGHJ1 {tb_blast.iloc[0].sacc} {sstart} {send} \n{jgene}\n")
