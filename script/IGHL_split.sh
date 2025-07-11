# extarct all seq 
# result/fullVDJ.fasta
# align them into 2 scaffolds
makeblastdb -parse_seqids -dbtype nucl -in data/2chrom.fa -out blastdb/2chrom

# blast all seq to 2 scaffolds  
blastn -query result/fullVDJ.fasta  -db blastdb/2chrom -evalue 1e-5 -max_target_seqs 1 -num_threads 60 -max_hsps 1 -outfmt '6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps' > result/blast2chrom.tsv 
