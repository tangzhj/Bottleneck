#coding=utf-8
import os
res_chrM="/public/home/chenbzh5/DB/hg19/human_g1k_v37_chrM/chrM.len"
res_ref="/public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/refdata-cellranger-atac-hg19-1.2.0/fasta/genome.fa"
merge_bam = "./pool_output/PBMC2_realign_NK/PBMC2_realign_NKcell.merge.bam"

file = open("./merge_NK2.sh", "w")
file.write("samtools merge " + merge_bam.replace(".bam", "A.bam") + " /public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/12-10_mk_realign_1000reads/PBMC_rep2/NK/*/*A.rmdup.bam &" + "\n")
file.write("samtools merge " + merge_bam.replace(".bam", "T.bam") + " /public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/12-10_mk_realign_1000reads/PBMC_rep2/NK/*/*T.rmdup.bam &" + "\n")
file.write("samtools merge " + merge_bam.replace(".bam", "G.bam") + " /public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/12-10_mk_realign_1000reads/PBMC_rep2/NK/*/*G.rmdup.bam &" + "\n")
file.write("samtools merge " + merge_bam.replace(".bam", "C.bam") + " /public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/12-10_mk_realign_1000reads/PBMC_rep2/NK/*/*C.rmdup.bam &" + "\n")
file.write("# samtools merge " + merge_bam + " ./pool_output/PBMC2_realign_NK/PBMC2_realign_NKcell*.bam" + "\n")
file.write("# samtools index " + merge_bam + "\n")
mpileup = merge_bam.replace(".bam", ".mt.pileup")
file.write("# samtools mpileup -l " + res_chrM + " -f " + res_ref + " -q 20 -Q 30 -x " + merge_bam + " > " + mpileup + "\n") 
snv = merge_bam.replace(".bam", ".snv")
file.write("# java -Xmx4g  -jar /public/home/chenbzh5/Tools/VarScan.v2.3.7.jar pileup2snp " + mpileup + " --min-var-freq 0.01  --min-reads2 2 > " + snv + "\n")
file.close()
os.system("sh merge_NK2.sh > merge_NK2.log 2>&1 &")
