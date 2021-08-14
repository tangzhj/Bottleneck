# Scripts of pipeline(Custom 10x Genomics scATAC-seq (mtscATAC-seq))

## 01_FastQC.sh
 Quality assessment for each "fq.gz"|"fastq.gz" files
## 02_Cellranger.sh
 Sequencing reads were aligned to a modified hg19 reference using CellRanger-ATAC count (v1.1.0)
## 03_Extract_chrM.sh
 Extract reads mapping in mitochondria DNA.
## 04_Split_bam.py
 Split bam file by cell barcord(tag'CR')
## 05_Single_Cell_SNV.py
 SNV calling for each cell.
## 06_Call_germline_mutation.py
 Bam files in each cell type were merged to call germline mutation (variant allele frequency > 90% was calculated as germline mutation in merged bam). 
## 07_Cell_SNP_Matrix.R
 Get matrix of mtDNA mutation. Rowname is cell barcode. Colname is variant site. 