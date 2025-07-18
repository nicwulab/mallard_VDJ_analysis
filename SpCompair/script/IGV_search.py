from Bio import SeqIO


fasta = 'fasta/IMGT000008/IMGT000008.fa'


for seq_record in SeqIO.parse(fasta, "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record.seq))

str(seq_record.seq[30975-1:31296].translate())
str(seq_record.seq[5347-1:5704+1])

str(seq_record.seq[5347-1:5704+1-9-12])



10365 - 10631
5347 - 5704

