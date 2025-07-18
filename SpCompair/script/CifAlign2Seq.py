from prody import parsePDB, matchAlign
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO



from prody import parsePDB, cealign




# 1) Load your mmCIF structures
prot1 = parsePDB('pdbModel/3jbw.cif')
prot2 = parsePDB('AF3/fold_rag1_rag2_human/fold_rag1_rag2_human_model_0.cif')

# 2) Select only Cα atoms
ca1 = prot1.select('name CA')
ca2 = prot2.select('name CA')


aln = cealign(pdb1, pdb2)



# 3) Perform structural alignment and get sequence alignment
pairs = matchAlign(ca1, ca2, seqid=0.3, overlap=0.5)

aln = pairs[0]
print(aln.getSeq1())

# 4) Extract the aligned sequences (with gaps) 
seq1 = aln.getSeq1()
seq2 = aln.getSeq2()

# 5) Create Biopython SeqRecord objects
rec1 = SeqRecord(Seq(seq1), id='prot1', description='Structure-based alignment')
rec2 = SeqRecord(Seq(seq2), id='prot2', description='Structure-based alignment')

# 6) Write the aligned sequences to a FASTA file
SeqIO.write([rec1, rec2], 'structure_based_alignment.fasta', 'fasta')
