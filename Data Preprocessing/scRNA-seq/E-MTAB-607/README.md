# Scripts of pipeline(Smart-seq2)

## 01_Single_Cell_SNV.py
 GATK pipline for mapping calling SNV(Input .fastq files, and output .snv files).
## 02_get_summary.sh
 Use pileup_inf_fj.pl to calculate average sequencing depth for each cell.
## 03_new_mtSNV_filter.sh
 Remove unconfident SNV(Every variant allele is 20 and more than 2 reads detected from both the forward and reverse orientation. retains mutations with balanced strand as (30% < forward / (forward +reverse)<70%).)
## 04_get_SNV_variation_linux
 Get variant allel frequence matrix of mtDNA mutation. Rowname is cell barcode. Colname is variant site.
## 05_get_SNV_count_linux.py
 Get reads count matrix of mtDNA mutation. Rowname is cell barcode. Colname is variant site.