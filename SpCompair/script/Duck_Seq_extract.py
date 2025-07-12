from Bio import SeqIO
from Bio.Seq import Seq

# read fasta
fasta_file = "/data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa"


with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        if record.id == 'ptg000068l':
            seq = str(record.seq)
            break

# Rag1
seq1 = seq[18171236:18174358]
seq2 = str(Seq(seq1).reverse_complement())
Rag1 = str(Seq(seq2).translate())

# Rag2
seq1 = seq[18160528:18162112]
Rag2 = str(Seq(seq1).translate())

with open("Rags.fa", "w") as f:
    f.write(f">Rag1\n{Rag1}\n")
    f.write(f">Rag2\n{Rag2}\n")






# Get Rag Sequences 
seq1 = seq[18171326:18173575]
seq2 = str(Seq(seq1).reverse_complement())
Rag1_ID = "ptg000189l|18171326:18173575|Rag-1|Duck"
Rag1 = str(Seq(seq2).translate())
Rag2_ID = "ptg000189l|18160528:18162112|Rag-2|Duck"
Rag2 = Rag2


# Get the RSS sequences

with open(fasta_file, "r") as file:
    for record in SeqIO.parse(file, "fasta"):
        if record.id == 'ptg000189l':
            seq = str(record.seq)
            break
        
# RSS for prediction
# IGLV

# 23 forward
F23_ID = "ptg000189l|534427:534472|Nicked 23-RSS intermediate forward strand|Duck"
F23 = seq[534427:534472]

R23_ID = "ptg000189l|534411:534472|Nicked 23-RSS intermediate reverse strand|Duck"
R23 = str(Seq(seq[534411:534472]).reverse_complement())

R12_ID = "ptg000189l|541038:541088|12-RSS signal end reverse strand|Duck" 
R12 = seq[541038: 541088]

F12_ID = "ptg000189l|541038:541072|12-RSS signal end forward strand|Duck"
F12 = str(Seq(seq[541038: 541072]).reverse_complement())

with open("AF3/Duck.fa", "w") as f:
    f.write(f">{Rag1_ID}\n{Rag1}\n")
    f.write(f">{Rag2_ID}\n{Rag2}\n")
    f.write(f">{F23_ID}\n{F23}\n")
    f.write(f">{R23_ID}\n{R23}\n")
    f.write(f">{R12_ID}\n{R12}\n")
    f.write(f">{F12_ID}\n{F12}\n")







