

# use biopython to get the seq NC_092617 from GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna

from Bio import SeqIO

for record in SeqIO.parse("data/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna", "fasta"):
    if record.id == "NC_092617.1":
        seq = record.seq
        print(seq[5874933:5875236])  # 0-based indexing
        break

print(seq[5874933:5875236+7+24+9])  # 0-based indexing
 
