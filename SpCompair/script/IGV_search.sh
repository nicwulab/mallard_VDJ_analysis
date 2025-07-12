# make a blast database from the fasta file

#!/bin/bash
makeblastdb -in fasta/IMGT000008/IMGT000008.fa \
  -out fasta/IMGT000008/IMGT000008 -dbtype nucl \
  -parse_seqids -title "IMGT000008"

tblastn -query fasta/chicken_IGLV1.pro \
  -db fasta/IMGT000008/IMGT000008 -outfmt 6 \
  -evalue 1e-2 | sort -n -k 9 > result/IMGT000008_IGLV1.tblastn

blastn -query data/IGLV_RSS.fa \
  -db fasta/IMGT000008/IMGT000008 -outfmt 6 \
  -evalue 10 -word_size 4| sort -n -k 9 > result/IMGT000008_RSS.blastn


cat  result/IMGT000008_IGLV1.tblastn  result/IMGT000008_RSS.blastn| sort -n -k 9 > result/IMGT000008_all.blast  


# find the RAG1 and RAG2 on Duck
tblastn -query fasta/Rag1_2ref.fa \
  -db /data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/Bird75_blast -outfmt 6 \
  -num_threads 60 > result/Duck_RAG1_2.tblastn
