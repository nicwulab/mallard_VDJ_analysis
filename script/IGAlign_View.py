#!/usr/bin/env python3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i','-I','--input')

args = parser.parse_args()
INPUT = args.input

with open(INPUT, 'r') as f:
    lines = f.read().split('\n')


Seqs = lines[1:]
Lines = len(Seqs)//8

print(lines[0])
for i in range(Lines):
    seq_a3 = Seqs[i*8]
    seq_a2 = Seqs[i*8+1]
    seq_a1 = Seqs[i*8+2]
    ref_aa = Seqs[i*8+3]
    ref_nt = Seqs[i*8+4]
    seq_nt = Seqs[i*8+5]
    seq_view = ""
    view_a1 = ""
    view_a2 = ""
    view_a3 = ""
    for ii in range(len(ref_nt)):
        if ref_nt[ii] != seq_nt[ii]:
            seq_view += f"\033[41;37m{seq_nt[ii]}\033[0m"
        else:
            seq_view += seq_nt[ii]
        if ref_aa[ii] != seq_a1[ii]:
            view_a1 += f"\033[41;37m{seq_a1[ii]}\033[0m"
        else:
            view_a1 += seq_a1[ii]
        a2 = seq_a2[ii]
        if seq_a2[ii] not in '- ' and len(ref_aa) < ii:
            if ref_aa[ii+1] == seq_a2[ii]:
                a2 = f"\033[41;37m{seq_a2[ii]}\033[0m"
            if ref_aa[ii-1] == seq_a2[ii]:
                a2 = f"\033[42;37m{seq_a2[ii]}\033[0m"
        view_a2 += a2
        a3 = seq_a3[ii]
        if seq_a3[ii] not in '- ' and len(ref_aa) < ii:
            if ref_aa[ii+1] == seq_a3[ii]:
                a3 = f"\033[41;37m{seq_a3[ii]}\033[0m"
            if ref_aa[ii-1] == seq_a3[ii]:
                a3 = f"\033[42;37m{seq_a3[ii]}\033[0m"
        view_a3 += a3
    print(view_a3)
    print(view_a2)
    print(view_a1)
    print(ref_aa)
    print(ref_nt)
    print(seq_view)
    print()




