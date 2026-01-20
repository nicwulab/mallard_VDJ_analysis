#!/bin/bash


# create the IGV/D/J fasta file for IgBlast database
cat  ../final/IG[HL]V.fa| awk '{print $1}' > fasta/IGV.fasta
cat  ../final/IG[HL]J.fa| awk '{print $1}' > fasta/IGJ.fasta
cat  ../final/IG[HL]D.fa| awk '{print $1}' > fasta/IGD.fasta

# create the IgBlast database
makeblastdb -parse_seqids -dbtype nucl -in fasta/IGV.fasta -out fasta/IGV
makeblastdb -parse_seqids -dbtype nucl -in fasta/IGJ.fasta -out fasta/IGJ
makeblastdb -parse_seqids -dbtype nucl -in fasta/IGD.fasta -out fasta/IGD

# create the aux file for IGJ
# the aux file is needed for IgBlast to recognize the J gene database
# the format is:
# <sequence_id>  start  Type  length
# The type is "JH" for heavy chain J genes and "JL" for light chain J genes

awk 'BEGIN{FS="\n";RS=">"}NR>1{split($1,a," ");seq_id=a[1];seq="";for(i=2;i<=NF;i++)seq=seq $i;len=length(seq);if(seq_id ~ /IGKJ/ || seq_id ~ /IGLJ/) type="JL"; else type="JH"; print ">"seq_id"\t1\t"type"\t"len}' fasta/IGJ.fasta| sed 's/>//;s/\t1\tJH/\t3\tJH/;s/\t1\tJL/\t2\tJL/' > fasta/Duck.aux

