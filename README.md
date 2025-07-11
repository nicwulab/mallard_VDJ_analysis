# Duck WGS VDJ annotation

location of data files:
`../data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa`

Index:
    - Bowtie2: `../data/20241202_DuckWGS_assemble/Bird75_bowtie`
    - Blast+: `../data/20241202_DuckWGS_assemble/Bird75_blast`

location of VDJ:
`../2411_4DcukHA/all_VDJ/*/outs/all_contig.fastq`
`../2411_4DcukHA/all_VDJ/*/outs/all_contig.fasta`


## Method:

1. Extract the seq from frame works and cdrs and align to the genome. After using cellranger to call VDJ, only complete sequences (has detected frame work region 1 to frame work region 4) would be used to trunked into the frame work sequence and cdr sequence. Then all of them are aligned to the genome to find the location of them by using blast.

2. summarizing the blast results for fr segments

## Plan

- [x] Extract the seq from frame works and cdrs and align to the genome
- [x] Iteration to find the V and J genes
- [ ] Identify the true VJ by signal sequences
- [ ] Identify the D gene by the junctions


So, After the general searching, it seems like the most os results are in the scaffold 'ptg000189l'. So, I will focus on this scaffold and run the blast again with less stringent parameters.


### align the tail of the V gene

```bash
awk -F"\t" '{print ">"$1"\n"$14}' result/potential_fr_complete.tsv> result/v_gene_tail.fa
muscle -align result/v_gene_tail.fa -output result/v_gene_tail.fasta
```






## VDJSeq-Solver: In Silico V(D)J Recombination Detection Tool

- mapping reads on rearranged gene segments poses a number of challenges:  V, D, J regions and by the inserted or deleted nucleotides due to enzymatic processes.

### Methods

VJ gene obtain:
    - bwotie to map the reads to the genome: unmmapepd reads are used to obtain the VJ gene segments
    - TopHat: to obtain the variable length alignments due to the presence of junction breakpoints (maybe I can try the TopHat-Fusion)




Some tools may helpful in this task:
- IMGT/V-QUEST: alignment based tool, germline sequences need to be provided (or in the database)
    1. Giudicelli V, Chaume D, Lefranc MP. IMGT/V-QUEST, an integrated software program for immunoglobulin and T cell receptor V-J and V-D-J rearrangement analysis. Nucleic Acids Research. 2004;32: W435–W440. pmid:15215425
    26. Brochet X, Lefranc MP, Giudicelli V. IMGT/V-QUEST: an algorithm for Immunoglobulin and T cell receptor sequence analysis. Actes des Journes Ouvertes Biologie, Informatique et Mathematiques. JOBIM 2007. 2007: 329–330.
    27. Giudicelli V, Brochet X, Lefranc MP. IMGT/V-QUEST: IMGT standardized analysis of the immunoglobulin (IG) and T cell receptor (TR) nucleotide sequences. Cold Spring Harb Protoc. 2011;6: pdb–prot5633. pmid:21632778
    28. Alamyar E, Duroux P, Lefranc MP, Giudicelli V. IMGT() tools for the nucleotide analysis of immunoglobulin (IG) and T cell receptor (TR) V-(D)-J repertoires, polymorphisms, and IG mutations: IMGT/V-QUEST and IMGT/HighV-QUEST for NGS. Methods Mol Biol. 2012;882: 569–604. pmid:22665256

- IMGT/JunctionAnalysis: The junction is here defined as the region starting at the second conserved cysteine of the V region at position 104 (2nd-CYS) and ending with the conserved tryptophan (J-TRP for the IGH chains) or the conserved phenylalanine (J-PHE for the IGL chains or the TCR chains) at position 118.
    *V*************[YF]
    1. Monod MY, Giudicelli V, Chaume D, Lefranc MP. IMGT/JunctionAnalysis: the first tool for the analysis of the immunoglobulin and T-cell receptor complex V-J and V-D-J JUNCTIONs. Bioinformatics. 2004;20: i379–i385.
    31. Giudicelli V, Lefranc MP. IMGT/junctionanalysis: IMGT standardized analysis of the V-J and V-D-J junctions of the rearranged immunoglobulins (IG) and T cell receptors (TR). Cold Spring Harb Protoc. 2011;6: 716–725.

- JOINSOLVER deals with the difficulty of D gene segment assignments giving a higher score for longer consecutive nucleotides matches, but searches for two relatively conserved motifs TAT TAC TGT and C TGG GG to find the extreme points of the Third Complementarity Determining Region (CDR3) that is the most variable part of BCR and TCR.
    1. 29.Souto-Carneriro MM, Longo NS, Russ DE, Sun HW, Lipsky PE. Characterization of the Human IG Heavy Chain Antigen Binding Complementarity Determining Region 3 Using a Newly Developed Software Algorithm, JOINSOLVER. The Journal of Immunology. 2004;172: 6790–6802.

- SoDA
    1. 32.Volpe JM, Cowell LG, Kepler TB. SoDA: implementation of a 3D alignment algorithm for inference of antigen receptor recombinations. Bioinformatics 2006;22: 438–444. pmid:16357034

-  VDJSolver [25] and iHMMune-align [24] apply instead statistical models to obtain the optimized parameters fitting to the rearranged sequence: The model robustness heavily depends on the quality and diversity of the training data sets
    1. 24.Gaeta BA, Malming HR, Jackson KJ, Bain ME, Wilson P, Collins AM. IHMMune-align: hidden Markov model-based alignment and Bioinformatics of germline genes in rearranged immunoglobulin gene sequences. Immunology. 2007;23: 1580–1587.
    25. Laursen O, Nielsen M, Larsen SR, Barington T. No evidence for the use of DIR, D-D fusions, chromosome 15 open reading frames or VH replacement in the peripheral repertoire was found on application of an improved algorithm, JointML, to 6329 human immunoglobulin H rearrangements. Immunology. 2006;119: 265–277.

