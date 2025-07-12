from Bio import SeqIO
from Bio.Seq import Seq

# read fasta
fasta_file = "fasta/IMGT000008/IMGT000008.fa"

with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        seq = str(record.seq)

# reverse complement
seq1 = seq[30975:31296]
seq2 = str(Seq(seq1).reverse_complement())

# RSS for prediction
# IGLV
seq[31257:31302]
str(Seq(seq[31257-16:31302]).reverse_complement())

# IGLJ
seq[33033:33083]
str(Seq(seq[33033:33083-16]).reverse_complement())
seq[33033:33083]





str(seq2)
