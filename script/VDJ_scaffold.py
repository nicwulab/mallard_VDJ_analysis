import os

fa_list = [i for i in os.popen('ls result/[fc]*.fa').read().split() if 'fr1'  in i][-1:]

for fa in fa_list:
    id = fa.split('/')[-1].split('.')[0]
    # checking if the result is already exist
    if os.path.exists(f"result/ptg000063l_{id}.tsv"):
        continue
    cmd = f"blastn -query {fa} -db ../data/20241202_DuckWGS_assemble/ptg000063l -evalue 1e-4 -max_target_seqs 100 -num_threads 60 -max_hsps 100 -outfmt '6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps' -word_size 15 -qcov_hsp_perc 100 -out result/ptg000063l_{id}.tsv"
    os.system(cmd)
    cmd = f"blastn -query {fa} -db ../data/20241202_DuckWGS_assemble/ptg000189l -evalue 1e-4 -max_target_seqs 100 -num_threads 60 -max_hsps 100 -outfmt '6 qacc sacc pident qcovs qlen qstart qend sstart send length mismatch gaps' -word_size 15 -qcov_hsp_perc 100 -out result/ptg000189l_{id}.tsv"
    os.system(cmd)
