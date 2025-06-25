# cut-run_repeat_analysis
Scripts to prepare data for DEseq2 analysis for cut&amp;run data mapping to centromeric repeat regions

This analysis accompanies the Sidhwani et al. (2025) paper "Histone H3 lysine methyltransferase activities control compartmentalization of human centromeres." Specifically, this pipeline is designed to take an adapter-trimmed cut&run dataset from 2x150 bp paired end Illumina sequencing and calculate statistically significant changes in reads mapping to various repeatmasker annotations within their centromeric contexts. This pipeline has been designed to handle "multimappers" i.e. repeat-derived reads that map to >1 location in the genome. For Straight lab users, all scripts and the bowtie1 dependancy are stored in `/home/groups/astraigh/shared/cutnrun_repeat_analysis/`

Please note, the python scripts used in this pipeline were written by Georgi K. Marinov and published in his original manuscript https://doi.org/10.1016/j.devcel.2015.01.013. They are also available on https://github.com/georgimarinov/GeorgiScripts

To run this pipeline, you will need:
1. adaptor trimmed fastq files from Illumina PE sequencing for a control and a treatment with a minimum of 2 replicates per condition
   NOTE: For ease of readibility, we will use the control examples WT_K9me3_r1, WT_K9me3_r2 and the treatment examples DKO_K9me3_r1 and DKO_K9me3_r2. All code where only WT_K9me3_r1 is mentioned must in fact be run for all 4 conditions/replicates
2. repeatmasker annotation bedfile (for most recent CHM13 version, see https://github.com/marbl/CHM13)
3. centromeric annotation bedfile (for most recent CHM13 version, see https://github.com/marbl/CHM13)
4. complete genome assembly fasta (preferably isogenomic assembly, but can use CHM13 from https://github.com/marbl/CHM13)
5. chromosome sizes txt file (https://github.com/marbl/CHM13)
6. Bash, R, python3

## Dependencies
You will need to install the following softwares:
1. BBTools (https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/)
2. bowtie1 patch that adds NH tags (https://mitra.stanford.edu/kundaje/marinovg/oak/programs/bowtie-1.0.1+hamrhein_nh_patch/)
3. samtools  (https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2)
4. bedtools (https://bedtools.readthedocs.io/en/latest/)
5. DEseq2 and its dependencies like tidyverse (https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

## Pipeline running instructions

### Step 0: Prepare fastq files
1. Use your favorite Illumina adaptor trimming method to trim sequencing adaptors from fastq files. We use Trimmomatic.
2. Prepare bowtie1 index for the T2Tv2.0 genome or another isogenomic sequence

### Step 1: Pairing reads with perfect overlap
With 2 x 150 bp PE sequencing, we first merged reads with no mismatch in the overlap. This is to get slightly longer fragments to align.
```
bbmerge.sh pfilter=1 in1={WT_K9me3_r1.PE1} in2={WT_K9me3_r1.PE2} out=WT_K9me3_r1.pfilter1.fastq
```
### Step 2: bowtie1 alignment and samtools index
1. Paired reads are aligned with either the strictest approach of no mismatches allowed and only unique alignments (-m 1 -v 0) or 1 mismatch allowed and all alignments reported (-v 1 -a). Of not, the --best and --strata options will choose for the best quality alignments even with -v 1 -a options. 

In our paper, DEseq2 analysis for heterochromatin marks H3K9me2, H3K9me3 and H3K27me3 was performed with -m 1 -v 0. However, using these options for epigenetic marks concentrated specifically at the centromere, like CENP-A, is not ideal. For CENP-A, we instead used -v 1 -a for further processing. Additionally, all visualization tracks for Integrated Genome Viewer (IGV), irrespective of the epigenetic mark, were made with the -v 1 -a options. We will only show steps for the -v1 -a files where they differ in processing from the uniquely aligned bams.

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
4. For chromosome-by-chromosome analysis to generated heatmaps, append the chr number of $1 to $4 of the file
   
```
awk -F'\t' 'BEGIN{OFS="\t"} {$4 = $4 "_" $1; print}' chm13v2.0_rmsk_censatanno.sorted.bed > chm13v2.0_rmsk_censatanno.chr.sorted.bed
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
We use custom python scripts (uploaded in scripts) for counting reads. Add the --uniqueBAM tag if the bam files were generated using -m 1 -v 0. If generated using -v 1 -a, do not add this tag. If bam files were generated using -v 1 -a, this script will normalize multimappers by the number of genomic loci they map to (based on the NH tag added by the patched bowtie1. So for example, if one read maps equally well to 10,000 loci, each locus is assigned a weighted score of 0.0001.

```
python bedRawReadCountsBAM.py3.py chm13v2.0_rmsk_censatanno.sorted.bed 0 WT_K9me3_r1.pfilter1.v0m1.bam WT_K9me3_r1.RM.censat.counts -uniqueBAM &&
```

```
python bedRawReadCountsBAM.py3.py hs1.rmsk_censat.complement.1kb.bed 0 WT_K9me3_r1.pfilter1.v0m1.bam WT_K9me3_r1.RM.complement.counts -uniqueBAM &&
```

### Step 6: Sum reads that fall within the same repeatmasker_censat category 
We use a custom python script (uploaded in scripts) for summing reads that fall within the same repeatmasker_censat category (e.g. ALR/Alpha_hor or, for chromosome-by-chromosome analysis, ALR/Alpha_hor_chr1). Make sure the column you use to base the summation on is column 4 (by bash) or column 3 (by python) 

```
python SumDups.py3.py WT_K9me3_r1.RM.censat.counts 3 7 WT_K9me3_r1.RM.censat.counts.sum
```

### Step 7: Merge counts from repeats to genomic complement and make a table 

```
cat <(awk 'BEGIN {OFS="\t"} {print "bin"NR, $4}' WT_K9me3_r1.RM.complement.counts)  <(awk 'BEGIN {OFS="\t"} {print $1, $2}' WT_K9me3_r1.censat.counts.sum) > WT_K9me3_r1.RM.censat.merge.counts.table
```

### Step 8: Make a DEseq2 counts table

```
sort -k1,1 WT_K9me3_r1.RM.censat.merge.counts.table | join -1 1 -2 1 -t$'\t' - <(sort -k1,1 WT_K9me3_r2.RM.censat.merge.counts.table) | join -1 1 -2 1 -t$'\t' - <(sort -k1,1 DKO_K9me3_r1.RM.censat.merge.counts.table) | join -1 1 -2 1 -t$'\t' - <(sort -k1,1 DKO_K9me3_r2.RM.censat.merge.counts.table) > K9me3_WT_r12_DKOr12.counts.table
```

### Step 9: If -v 1 -a options were used, integerize counts

```
awk 'BEGIN{FS=OFS="\t"} {print $1, sprintf("%.0f", $2), sprintf("%.0f", $3), sprintf("%.0f", $4), sprintf("%.0f", $5)}' K9me3_WT_r12_DKOr12.counts.table > K9me3_WT_r12_DKOr12.counts.int.table
```

### Step 10: Add a header to the table

```
echo -e "\tWT_r1\tWT_r2\tDKO_r1\tDKO_r2" | cat - K9me3_WT_r12_DKOr12.counts.table > K9me3_WT_r12_DKOr12.counts.htable
```
### Step 11: Make a sample metadata file
1. We use vi to make the metadata file
```
vi K9me3_WT_r12_DKOr12.meta
```
2. Insert the following, making sure the names match with the table header from Step 10.

```
	sampletype
WT_r1	control
WT_r2	control
DKO_r1	knockout
DKO_r2	knockout
```

### Step 12: Make a DEseq2 script
1. We use vi to make the DEseq2 script
```
vi DEseq2.R
```
2. Insert the following

```
library(DESeq2)
library(tidyverse)

K9me3_m1v0_censat <- read.table("K9me3_WT_r12_DKOr12.counts.htable", header = TRUE, sep = "\t", row.names = 1)

metadata <- read.table("K9me3_WT_r12_DKOr12.meta", header = TRUE, sep = "\t", row.names = 1)

print(ncol(K9me3_m1v0_censat))
print(nrow(metadata))

print(rownames(metadata))
print(colnames(K9me3_m1v0_censat))

dds <- DESeqDataSetFromMatrix(K9me3_m1v0_censat, metadata, design = ~ sampletype)

dds <- DESeq(dds)

contrast_K9me3_m1v0_censat <- c("sampletype", "knockout", "control")

res_table <- results(dds, contrast=contrast_K9me3_m1v0_censat, alpha = 0.01)

res_table <- lfcShrink(dds, contrast=contrast_K9me3_m1v0_censat, type="ashr", res=res_table)

summary(res_table, alpha = 0.01)

padj.cutoff <- 0.01

res_table_tb <- res_table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig <- res_table_tb %>%
        filter(padj < padj.cutoff)
saveRDS(sig, file = "sig_WTr12_DKOr12_K9me3_m1v0_censat.rds")
saveRDS(res_table_tb, file = "res_table_tb_WTr12_DKO_r12_K9me3_m1v0_censat.rds")
write.table(sig, file = "K9me3_censat_WTr12_DKOr12_sig.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

And voila! You have a table containing your significant hits!

### Step 13: Creating bigwig tracks for visualization
Use a custom python script (in scripts) to generate normalized bedGraph tracks for visualization. This script takes a bam file containing NH tags and generated with -v 1 -a to output a bedGraph where reads are nomalized to the number of genomic loci they map to.

```
python bedRPKMfrombam.py3.py hs1.1kb.windows_lastbinremoved.bed 0 WT_K9me3_r1.pfilter1.v1a.bam hs1.chromsizes.txt WT_K9me3_r1.pfilter1.v1a.bdg -excludeReadsMappingToOtherChromosomes -RPM
```










