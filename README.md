# cut-run_repeat_analysis
Scripts to prepare data for DEseq2 analysis for cut&amp;run data mapping to centromeric repeat regions

This analysis accompanies the Sidhwani et al. (2025) paper "Histone H3 lysine methyltransferase activities control compartmentalization of human centromeres." Specifically, this pipeline is designed to take an adapter-trimmed cut&run dataset from 2x150 bp paired end Illumina sequencing and calculate statistically significant changes in reads mapping to various repeatmasker annotations within their centromeric contexts. This pipeline has been designed to handle "multimappers" i.e. repeat-derived reads that map to >1 location in the genome. For Straight lab users, all scripts and the bowtie1 dependancy are stored in '/home/groups/astraigh/shared/cutnrun_repeat_analysis/'

To run this pipeline, you will need:
1. adaptor trimmed fastq files from Illumina PE sequencing for a control and a treatment with a minimum of 2 replicates per condition
2. repeatmasker annotation bedfile
3. centromeric annotation bedfile
4. complete genome assembly (preferably isogenomic assembly
5. Bash, R, python3

## Dependencies
You will need to install the following softwares:
1. bowtie1 patch that adds NH tags (https://mitra.stanford.edu/kundaje/marinovg/oak/programs/bowtie-1.0.1+hamrhein_nh_patch/)
2. samtools  (https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2)
3. bedtools (https://bedtools.readthedocs.io/en/latest/)
4. DEseq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)



