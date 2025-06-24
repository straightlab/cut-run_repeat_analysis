# cut-run_repeat_analysis
Scripts to prepare data for DEseq2 analysis for cut&amp;run data mapping to centromeric repeat regions

This analysis accompanies the Sidhwani et al. (2025) paper "Histone H3 lysine methyltransferase activities control compartmentalization of human centromeres." Specifically, this pipeline is designed to take an adapter-trimmed cut&run dataset from 2x150 bp paired end Illumina sequencing and calculate statistically significant changes in reads mapping to various repeatmasker annotations within their centromeric contexts. This pipeline has been designed to handle "multimappers" i.e. repeat-derived reads that map to >1 location in the genome. For Straight lab users, all scripts and the bowtie1 dependancy are stored in `/home/groups/astraigh/shared/cutnrun_repeat_analysis/`

To run this pipeline, you will need:
1. adaptor trimmed fastq files from Illumina PE sequencing for a control and a treatment with a minimum of 2 replicates per condition
   NOTE: For ease of readibility, we will use the control examples WT_K9me3_r1, WT_K9me3_r2 and the treatment examples DKO_K9me3_r1 and DKO_K9me3_r2. All code where only WT_K9me3_r1 is mentioned must in fact be run for all 4 conditions/replicates
2. repeatmasker annotation bedfile
3. centromeric annotation bedfile
4. complete genome assembly fasta (preferably isogenomic assembly)
5. Bash, R, python3

## Dependencies
You will need to install the following softwares:
1. BBTools (https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/)
2. bowtie1 patch that adds NH tags (https://mitra.stanford.edu/kundaje/marinovg/oak/programs/bowtie-1.0.1+hamrhein_nh_patch/)
3. samtools  (https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2)
4. bedtools (https://bedtools.readthedocs.io/en/latest/)
5. DEseq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

## Pipeline running instructions

### Step 0: Prepare fastq files
1. Use your favorite Illumina adaptor trimming method to trim sequencing adaptors from fastq files. We use Trimmomatic.
2. Prepare bowtie1 index for the T2Tv2.0 genome

### Step 1: Pairing reads with perfect overlap
With 2 x 150 bp PE sequencing, we first merged reads with no mismatch in the overlap. This is to get slightly longer fragments to align.
```
bbmerge.sh pfilter=1 in1={WT_K9me3_r1.PE1} in2={WT_K9me3_r1.PE2} out=WT_K9me3_r1.pfilter1.fastq
```
### Step 2: bowtie1 alignment and samtools index
1. Paired reads are aligned with either the strictest approach of no mismatches allowed and only unique alignments (-m 1 -v 0) or 1 mismatch allowed and all alignments reported (-v 1 -a). Of not, the --best and --strata options will choose for the best quality alignments even with -v 1 -a options. 

In our paper, DEseq2 analysis for heterochromatin marks H3K9me2, H3K9me3 and H3K27me3 was performed with -m 1 -v 0. However, using these options for epigenetic marks concentrated specifically at the centromere, like CENP-A, is not ideal. For CENP-A, we instead used -v 1 -a for further processing. Additionally, all IGV plots, irrespective of the epigenetic mark, were made with the -v 1 -a options. We will only show steps for the -v1 -a files where they differ in processing from the uniquely aligned bams.

```
bowtie <path to bowtie index> -p 12 -t -v 0 -m 1 --best --strata -q --sam-nh -X 1000 --sam WT_K9me3_r1.pfilter1.fastq -S | samtools view -bS | samtools sort > WT_K9me3_r1.v0m1.bam
```
OR

```
bowtie <path to bowtie index> -p 12 -t -v 1 -a --best --strata -q --sam-nh -X 1000 --sam WT_K9me3_r1.pfilter1.fastq -S | samtools view -bS | samtools sort > WT_K9me3_r1.v1a.bam
```
2. indexing bam files
   
```
samtools index -@ 16 WT_K9me3_r1.v0m1.bam
```

### Step 3: Prepare repeat_censat bedfile for read counting
In our paper, we used a bedfile where the repeatmasker annotations were intersected with centromeric annotations. 
1. To prepare this bedfile, we first prepare a bedfile that contains information from both the repeatmasker file and the centromeric annotation file

```
bedtools intersect -a repeatmasker.bed -b centsatannotation.bed > repeatmasker_censat.bed -wa -wb
```
2. Next, we we used append the centromeric annotation from column 13 to the repeatmasker annotation in column 4 and subsequently kept only the first 7 columns. In this way, column 4 of the new bed file contains information about both the repeatmasker and the centromeric annotations for the genomic coordinates in columns 1-3.

```
awk -F'\t' 'BEGIN{OFS="\t"} {
    split($13, a, "_");
    print $1, $2, $3, $4 "_" a[1], $5, $6, $7
}' repeatmasker_censat.bed > chm13v2.0_rmsk_censatanno.bed
```
3. Sort bedfile if needed

```
bedtools sort -i chm13v2.0_rmsk_censatanno.bed > chm13v2.0_rmsk_censatanno.sorted.bed
```

### Step 4: Prepare genome bedfile for read counting
In addition to counting reads aligning to repeat regions, we additionally need to count reads aligning to genomic bins. We choose 1 kb genomic bins, however the size can be changed depending on your application.
1. First, get the genomic complement of the repeatmasker_censat file
```
bedtools complement -i chm13v2.0_rmsk_censatanno.sorted.bed -g genome.fa > hs1.rmsk_censat.complement.bed
```
2. Then, make 1 kb windows

```
bedtools makewindows -b hs1.rmsk_censat.complement.bed -w 1000 > hs1.rmsk_censat.complement.1kb.bed
```
### Step 5: Count reads in the repeatmasker_censat regions and the genomic complement
We use custom python scripts (uploaded in scripts) for counting reads. Add the --uniqueBAM tag if the bam files were generated using -m 1 -v 0. If generated using -v 1 -a, do not add this tag.

```
python bedRawReadCountsBAM.py3.py chm13v2.0_rmsk_censatanno.sorted.bed 0 WT_K9me3_r1.pfilter1.v0m1.bam WT_K9me3_r1.RM.censat.counts -uniqueBAM &&
```

```
python bedRawReadCountsBAM.py3.py hs1.rmsk_censat.complement.1kb.bed 0 WT_K9me3_r1.pfilter1.v0m1.bam WT_K9me3_r1.RM.censat.complement.counts -uniqueBAM &&
```








