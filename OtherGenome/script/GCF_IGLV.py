# read sequence 'NC_092602' from data/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna
from Bio import SeqIO
from Bio.Seq import Seq

def read_sequence(file_path, sequence_id):
    for record in SeqIO.parse(file_path, "fasta"):
        if record.id == sequence_id:
            return str(record.seq)
    return None

# the column is qacc sacc pident qcovs qlen qstart qend sstart send
def read_blast_results(file_path, threshold=95):
    blast_hits = []
    with open(file_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            hit = {
                "qacc": parts[0],
                "sacc": parts[1],
                "pident": float(parts[2]),
                "qcovs": float(parts[3]),
                "qlen": int(parts[4]),
                "qstart": int(parts[5]),
                "qend": int(parts[6]),
                "sstart": int(parts[7]),
                "send": int(parts[8]),
            }
            blast_hits.append(hit)
    # filter the hits with  qcovs >= 95
    filtered_hits = [hit for hit in blast_hits if hit["qcovs"] >= threshold]
    return filtered_hits

def ScoreCalc(fw):
    score = 0
    for i in range(len("CACAGTG")):
        if fw[i] == "CACAGTG"[i]:
            score += 1
    for i in range(len("CAAAAACC")):
        if fw[-8:][i] == "CAAAAACC"[i]:
            score += 1
    return score


# Get the Sequence 
file_path = "data/GCF_047663525.1_IASCAAS_PekinDuck_T2T_genomic.fna"
sequence_id = "NC_092602.1"
sequence = read_sequence(file_path, sequence_id)


# read the blast results 'result/GCF_IGLJ.blast.tsv'

blast_file = "result/GCF_IGLJ.blast.tsv"
filtered_hits = read_blast_results(blast_file) 

# get teh IGLJ
IGVJ = filtered_hits[0]
IGVJ_seq = sequence[IGVJ["sstart"]-29:IGVJ["send"]]
print("IGHJ sequence:", IGVJ_seq)

# Get all potential IGLV
blast_file_IGLV = "result/GCF_IGLV.blast.tsv"
filtered_hits_IGLV = read_blast_results(blast_file_IGLV)
chrom = "NC_092602"
filtered_hits_IGLV_chrom = [hit for hit in filtered_hits_IGLV if hit["sacc"] == chrom]

Seq_fw = []
Seq_rv = []
for i, hit in enumerate(filtered_hits_IGLV_chrom):
    if hit["send"] > hit["sstart"]:
        for steps in range(-10,10):
            fw = sequence[hit["send"]+1+steps:hit["send"]+40+steps].upper()
            score = ScoreCalc(fw)
            if score >= 12:
                Seq_fw.append(IGLV_seq)
                print(f"IGLV sequence {i+1} (step {steps}):", fw, score)
        if "CACAGTG" in IGLV_seq:
            if "CAAAAACC" in IGLV_seq:
                print(f"IGLV sequence {i+1}:", IGLV_seq)
    else:
        for steps in range(-10,10):
            rv = sequence[hit["sstart"]-40+steps:hit["sstart"]-1+steps].upper()
            rv = str(Seq(rv).reverse_complement())
            score = ScoreCalc(rv)
            if score >= 12:
                Seq_rv.append(rv)
                print(f"IGLV sequence {i+1} (step {steps}):", rv, score, "rev")

        IGLV_seq = sequence[hit["sstart"]-40:hit["sstart"]-1].upper()
        IGLV_seq = str(Seq(IGLV_seq).reverse_complement())
        Seq_rv.append(IGLV_seq)

Seq_fw_dict = {}

ScoreCalc(fw)

for fw in Seq_fw:
    score = 0
    for i in range(len("CACAGTG")):
        if fw[i] == "CACAGTG"[i]:
            score += 1
    for i in range(len("CAAAAACC")):
        if fw[-8:][i] == "CAAAAACC"[i]:
            score += 1
    Seq_fw_dict[fw] = score


for fw,score in Seq_fw_dict.items():
    if score >= 8:
        print(fw, score)


# write the seq as txt, no name need
with open("result/IGLV_forward.txt", "w") as f:
    for seq in Seq_fw:
        f.write(seq + "\n")






