# Methods: IGH Pseudogene Complete Annotation Pipeline

## Overview

This pipeline extracts and annotates complete IGH (Immunoglobulin Heavy chain) pseudogene sequences from the duck genome assembly. The method addresses limitations of previous approaches by:
1. Using functional gene sequences to define complete pseudogene heads
2. Using D gene and V gene alignments to define pseudogene tails
3. Performing comprehensive sequence alignment and validation

## Pipeline Workflow

### Step 1: Data Preparation

#### 1.1 Load Genome Assembly
- **Input**: Pre-prepared genome assembly file
  - Path: `/data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa`
  - Format: FASTA format containing assembled contigs
  - **Status**: Pre-prepared (from genome assembly pipeline)

#### 1.2 Extract VDJ Sequences
- **Input**: VDJ annotation files from previous analysis
  - Source: `../2411_4DcukHA/result/*VDJ.tsv`
  - **Status**: Generated from previous VDJ annotation pipeline (CellRanger or similar)
- **Process**: Extract full VDJ sequences from TSV files and convert to FASTA format
- **Output**: `result/fullVDJ.fasta` (all VDJ sequences)

#### 1.3 Prepare Scaffold Database
- **Input**: Two-chromosome scaffold file
  - Path: `data/2chrom.fa`
  - Contains: `ptg000189l` (IGH scaffold) and `ptg000063l` (IGL scaffold)
  - **Status**: Pre-prepared (extracted from genome assembly)
- **Process**: Create BLAST database for the two scaffolds
- **Output**: `blastdb/2chrom` (BLAST database)

#### 1.4 Map VDJ Sequences to Scaffolds
- **Input**: `result/fullVDJ.fasta`
- **Process**: BLASTN alignment of all VDJ sequences to the two-chromosome database
  - Parameters: `evalue 1e-5`, `max_target_seqs 1`, `max_hsps 1`
- **Output**: `result/blast2chrom.tsv` (BLAST results)

#### 1.5 Separate IGH and IGL Sequences
- **Input**: `result/blast2chrom.tsv`
- **Process**: 
  - Filter sequences mapping to `ptg000189l` → IGH sequences
  - Filter sequences mapping to `ptg000063l` → IGL sequences
- **Output**: 
  - `result/IGH.list` (IGH sequence IDs)
  - `result/IGL.list` (IGL sequence IDs)
  - `result/IGH.fa` (IGH sequences in FASTA)
  - `result/IGL.fa` (IGL sequences in FASTA)

### Step 2: D Gene Extraction

#### 2.1 Extract D Gene Sequences
- **Input**: Genome assembly (`seq_dict`) and predefined D gene positions
- **Known D Gene Positions** (relative positions on scaffold):
  - IGHD1: (2083, 2160)
  - IGHD2: (3162, 3244)
  - IGHD3: (4132, 4213)
  - IGHD4: (5282, 5366)
  - IGHD5: (6190, 6281)
- **Known D Gene Sequences** (for validation):
  - IGHD1: `AGTGCTGGTGCTGGTTGGGTTA`
  - IGHD2: `GGTACTGGTGGTTGGGGTGCTAATGCA`
  - IGHD3: `AGTACTGGTGCTTGTTGTGGTGGTTA`
  - IGHD4: `GGTGCTGGTAGTGGTTACGGTGCTTATTA`
  - IGHD5: `GGTTATGCTAGTTGTGGTGGTTATACTTGTGCTTAT`
- **Process**: Extract D gene sequences from genome using known positions with offset adjustments
- **Output**: `final/IGHD.fa` (D gene sequences)

#### 2.2 Annotate D Genes in VDJ Sequences
- **Input**: 
  - `result/IGH.fa` (IGH VDJ sequences)
  - `final/IGHD.fa` (D gene sequences)
- **Process**: 
  - Create BLAST database for D genes
  - BLASTN alignment of IGH sequences to D genes
  - Parameters: `strand plus`, `evalue 1`, `word_size 4`
- **Output**: `result/IGH-IGHD.tsv` (D gene alignment results)

#### 2.3 Extract Leader Sequences (V Gene Heads)
- **Input**: `result/IGH-IGHD.tsv` (D gene alignment results)
- **Process**: Extract sequences upstream of D gene alignment start positions
  - For each VDJ sequence, extract the portion before the D gene match
  - This represents the leader/V gene head region
- **Output**: `result/Leader_IGHV.fa` (potential V gene leader sequences)

### Step 3: V Gene Extraction and Completion

#### 3.1 Functional V Gene Reference
- **Input**: Known functional IGHV1 sequence (first 99 nucleotides)
  - Sequence: `GCTGCCACCTTGGATGAGTCCGGAGGGGGCCTCGTGAGTCCCGGGGGGTCCCTGACCCTGGTCTGCAAGGGCTCCGGATACACCTTCAGCAGCTACGGCATGGGATGGGTGCGACAGGCACCCGGGAAGGGGCTCGAGTACGTCGCGAGTATTAACAGCAGTGGTAGTAGCACTTACTACGCGCCGGCGGTGAAGGGACGCTTCACCATCTCCAGGAACAACGGGCAGAGCACGCTCACCCTGCAGATGAACAGCCTCAAGGCCGAAGACACCGCCACCTACTACTGCGCGAAAGCTGCTGGT`
- **Process**: Create reference sequence file with first 99 nucleotides
- **Output**: `result/IGHV1-1_99.fa`

#### 3.2 Identify V Gene Heads
- **Input**: 
  - `result/IGHV1-1_99.fa` (functional V gene reference)
  - `result/Leader_IGHV.fa` (leader sequences)
- **Process**: 
  - Create BLAST database for leader sequences
  - BLASTN alignment of functional V gene to leader sequences
  - Parameters: `strand plus`, `evalue 1`, `qcov_hsp_perc 70`, `max_target_seqs 100000`
- **Output**: `result/IGHV1-1_99-Leader_IGHV.tsv` (alignment results)

#### 3.3 Extract Complete V Gene Sequences
- **Input**: `result/IGHV1-1_99-Leader_IGHV.tsv` (alignment results)
- **Process**: 
  - For each alignment, extract the complete V gene sequence
  - Start position: `sstart - qstart` (align leader start with functional V start)
  - Filter sequences with length ≥ 200 nucleotides
- **Output**: `result/IGHV.fa` (complete V gene sequences)

#### 3.4 Map V Genes to Genome Scaffold
- **Input**: 
  - `result/IGHV.fa` (V gene sequences)
  - Genome scaffold: `ptg000189l`
- **Process**: 
  - BLASTN alignment of V genes back to the IGH scaffold
  - Parameters: `evalue 1e-5`, `max_target_seqs 100`, `qcov_hsp_perc 70`
  - Include subject sequence in output (`sseq` field)
- **Output**: `result/IGHV-ptg000189l.tsv` (V gene to scaffold alignments)

#### 3.5 Refine V Gene Boundaries
- **Input**: `result/IGHV-ptg000189l.tsv` (alignment results)
- **Process**: 
  - Group alignments by genomic position
  - Identify distinct V gene clusters using position-based sorting
  - For each cluster:
    - Select alignments with `qcovs >= 100`
    - Determine MIN and MAX genomic coordinates
    - Perform protein-level alignment to find optimal reading frame
    - Test all 6 reading frames (3 forward, 3 reverse)
    - Use BLOSUM62 substitution matrix for alignment scoring
    - Select best alignment based on score
    - Adjust boundaries to remove gaps at sequence head
    - Filter sequences with alignment score > 50
- **Output**: 
  - `final/IGHV.fa` (nucleotide sequences)
  - `final/IGHV.pro` (protein sequences)
  - `result/IGHV_align/IGHV1-*.aln` (individual alignment files)

#### 3.6 Multiple Sequence Alignment
- **Input**: 
  - `final/IGHV.fa` (nucleotide sequences)
  - `final/IGHV.pro` (protein sequences)
- **Process**: MUSCLE alignment
- **Output**: 
  - `result/IGHV.aln` (nucleotide alignment)
  - `result/IGHV.pro.aln` (protein alignment)

### Step 4: J Gene Extraction

#### 4.1 Extract J Gene Sequence
- **Input**: Genome assembly and known J gene position
- **Known J Gene Position**: 
  - Start: 541834
  - End: 541928
- **Process**: Extract J gene sequence from genome using known position with offset adjustments
- **Output**: `final/IGHJ.fa` (J gene sequence)

## Input Files Summary

### Pre-prepared Files (External Dependencies)

1. **Genome Assembly**
   - Path: `/data3/wenkanl2/Tomas/data/20241202_DuckWGS_assemble/Bird75_min1k_trimmed_l0_cov90.p_ctg.fa`
   - Source: Genome assembly pipeline
   - Description: Complete genome assembly in FASTA format

2. **VDJ Annotation Files**
   - Path: `../2411_4DcukHA/result/*VDJ.tsv`
   - Source: Previous VDJ annotation analysis (likely CellRanger)
   - Description: TSV files containing VDJ sequence annotations with framework regions and CDRs

3. **Two-Chromosome Scaffold File**
   - Path: `data/2chrom.fa`
   - Source: Extracted from genome assembly
   - Description: Contains IGH scaffold (`ptg000189l`) and IGL scaffold (`ptg000063l`)

4. **IGH Scaffold File**
   - Path: `../data/20241202_DuckWGS_assemble/ptg000189l` or `data/ptg000189l.fa`
   - Source: Extracted from genome assembly
   - Description: IGH locus scaffold for detailed V gene mapping

### Generated Files (Intermediate)

1. **VDJ Sequences**: `result/fullVDJ.fasta` (from Step 1.2)
2. **BLAST Results**: `result/blast2chrom.tsv` (from Step 1.4)
3. **IGH Sequences**: `result/IGH.fa` (from Step 1.5)
4. **D Gene Sequences**: `final/IGHD.fa` (from Step 2.1)
5. **D Gene Alignments**: `result/IGH-IGHD.tsv` (from Step 2.2)
6. **Leader Sequences**: `result/Leader_IGHV.fa` (from Step 2.3)
7. **V Gene Alignments**: `result/IGHV1-1_99-Leader_IGHV.tsv` (from Step 3.2)
8. **V Gene Sequences**: `result/IGHV.fa` (from Step 3.3)
9. **V Gene to Scaffold Alignments**: `result/IGHV-ptg000189l.tsv` (from Step 3.4)

### Final Output Files

1. **IGHV Nucleotide Sequences**: `final/IGHV.fa`
2. **IGHV Protein Sequences**: `final/IGHV.pro`
3. **IGHV Alignments**: `result/IGHV.aln` (nucleotide), `result/IGHV.pro.aln` (protein)
4. **Individual V Gene Alignments**: `result/IGHV_align/IGHV1-*.aln`
5. **IGHD Sequences**: `final/IGHD.fa`
6. **IGHJ Sequence**: `final/IGHJ.fa`

## Scientific Methods

### 1. Pseudogene Head Completion

**Problem**: Previous methods failed to capture complete pseudogene heads, resulting in truncated sequences.

**Solution**: 
- Use a known functional V gene sequence (IGHV1) as a reference
- Extract leader sequences from VDJ sequences by identifying D gene positions
- Align functional V gene to leader sequences to identify the start position
- Extract complete V gene sequences starting from the aligned position

**Algorithm**:
1. For each VDJ sequence, identify D gene position using BLAST alignment
2. Extract sequence upstream of D gene as potential leader/V head
3. Align functional V gene (first 99 nt) to leader sequences
4. If alignment found (`qcov >= 70%`), extract complete V gene from position `sstart - qstart`

### 2. Pseudogene Tail Definition

**Problem**: Previous methods did not properly define pseudogene tails.

**Solution**:
- Map extracted V gene sequences back to the genome scaffold
- Use protein-level alignment to identify optimal reading frames
- Refine boundaries by removing gaps and adjusting for frame shifts

**Algorithm**:
1. BLAST V gene sequences to IGH scaffold (`ptg000189l`)
2. Group alignments by genomic position
3. For each position cluster:
   - Test all 6 reading frames (3 forward, 3 reverse)
   - Perform protein-level alignment using BLOSUM62 matrix
   - Select frame with highest alignment score
   - Adjust boundaries to remove leading gaps
   - Filter by alignment score threshold (>50)

### 3. Reading Frame Detection

**Method**: Multi-frame protein alignment
- Translate genomic sequence in all 6 possible reading frames
- Align translated sequences to functional V gene protein sequence
- Use BLOSUM62 substitution matrix for scoring
- Select frame with maximum alignment score
- Adjust genomic coordinates based on selected frame and gap positions

**Parameters**:
- Gap open penalty: -10
- Gap extension penalty: -0.5
- Substitution matrix: BLOSUM62

### 4. D Gene Identification

**Method**: Position-based extraction with known sequences
- Use predefined genomic positions for D genes (identified in previous analysis)
- Extract sequences with offset adjustments for signal sequences
- Validate extracted sequences against known D gene sequences

### 5. Sequence Alignment and Validation

**Multiple Sequence Alignment**:
- Use MUSCLE for both nucleotide and protein alignments
- Generate alignment files for downstream analysis and visualization

**Quality Control**:
- Filter V genes by minimum length (≥200 nucleotides)
- Filter alignments by query coverage (≥70% or ≥100% depending on step)
- Filter by alignment score threshold (>50)

### 6. Coordinate System

**Genomic Coordinates**:
- All positions are relative to the IGH scaffold (`ptg000189l`)
- D gene positions include offset adjustments: `position + 534855 + adjustments`
- J gene position: 541834-541928 (with offset adjustments)

## Technical Notes

### Dependencies
- **Python packages**: BioPython, pandas, numpy
- **External tools**: BLAST+ (makeblastdb, blastn), MUSCLE, seqkit
- **Data formats**: FASTA, TSV

### Known Issues
- Variables `fun_gene` and `tb_blast` are referenced but not defined in the script
  - These likely need to be set from previous analysis results
  - `fun_gene.sacc` should be `'ptg000189l'`
  - `tb_blast` should be a pandas DataFrame with scaffold information

### Pipeline Execution Order
1. Run data preparation steps (1.1-1.5)
2. Extract D genes (2.1-2.3)
3. Extract and refine V genes (3.1-3.6)
4. Extract J genes (4.1)

The pipeline can be run sequentially, with each step depending on outputs from previous steps.

---

## Methods (Condensed Description)

All IGH pseudogene annotation steps were performed in an integrated workflow to minimize fragmentation and redundancy. The duck genome assembly (Bird75_min1k_trimmed_l0_cov90.p_ctg.fa) served as the reference. VDJ sequences were extracted from four duck single-cell VDJ-Seq (scVDJ-Seq) assembly TSV files, converted to FASTA format (63144 in total), and mapped to the genome using BLASTN (e-value ≤ 1×10⁻⁵, max_target_seqs=1). IGH- and IGL-mapping sequences were identified by BLAST hit to the respective scaffolds (ptg000189l for IGH, ptg000063l for IGL).

D gene segments (IGHD1-5) were defined by manual identification of RSS (recombination signal sequence) motifs on the IGH scaffold, using imperfect but targeted searches to locate RSS-flanked segments. The resulting scaffold coordinates (2083-2160, 3162-3244, 4132-4213, 5282-5366, 6190-6281 relative to position 534855) reflect these empirically determined boundaries. Extracted D sequences included additional offset nucleotides to ensure complete RSS coverage and were validated against known duck D gene references. A BLASTN database of D gene sequences was constructed to profile D gene usage in VDJ repertoires and show that 8159 of 9255 IGH can be annotated. **88.16%** D gene were found.

To recover full V gene heads, the region upstream of each D gene alignment in VDJ sequences was considered the leader. Each leader was aligned to a functional IGHV1 gene reference (99 nt) with BLASTN (e-value=1, query coverage ≥70%). For each leader hit, a putative full-length V gene sequence was extracted, extending from the inferred start. Sequences <200 nt were excluded.

V gene pseudogene boundaries and reading frame were refined by mapping candidate V genes back to the IGH scaffold (BLASTN: e-value ≤ 1×10⁻⁵, query coverage ≥70%). Hits were grouped by genomic location. For each group, alignments with query coverage ≥100% were used to define cluster coordinates. Cluster-spanning genomic DNA was subjected to six-frame translation, with each frame aligned to the functional IGHV1 protein (pairwise local alignment, BLOSUM62, gap open -10, gap extension -0.5). The frame with the highest score was taken as the correct reading frame; leading gaps/frameshifts were trimmed to maximize coding sequence. Sequences with alignment score ≤50 were discarded.

The J gene was similarly extracted from a fixed IGH coordinate (541834-541928, with RSS adjustment) using the same extraction method as for D genes.

Final V and D gene sequences, both nucleotide and protein, were aligned using MUSCLE for further analysis. At all stages, sequences failing minimum criteria (length ≥200 nt, alignment coverage ≥70–100%, protein alignment score >50, frame consistency) were excluded. Whole-pipeline sequence and analysis management used BLAST+ (v2.x), MUSCLE, BioPython (for BLOSUM62 and translation), pandas, numpy, and seqkit as needed.

Pseudogene annotations required all these quality and statistical criteria: full-length segment recovery (by coverage), protein-level match to the functional gene, adequate sequence length, and frame integrity. Only sequences passing all criteria were kept for subsequent annotation and analysis.


