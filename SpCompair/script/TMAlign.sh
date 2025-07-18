
for tg in `ls AF3/*/*.cif`; do
  NAME=$(echo $tg| awk -F"/" '{print $NF}'| sed 's/\.cif/.pdb/g')
  TMalign pdbModel/3jbw.cif $tg -o result/TMalign/$NAME -outfmt 1 > result/TMalign/$NAME.fasta
done

