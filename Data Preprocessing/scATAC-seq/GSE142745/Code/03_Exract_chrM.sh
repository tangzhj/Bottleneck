#PBS -N scATAC
#PBS -j oe
#PBS -q batch
#PBS -o example.stdout
#PBS -e example2.stdout
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=2
#PBS -l mem=10000m
samtools view -b possorted_bam.bam chrM  > possorted_bam.chrM.bam
samtools view -h possorted_bam.chrM.bam > 1.sam
sed -i "s/chrM/chrMT/g"  1.sam
sed -i "27,87d"  1.sam
samtools view -bS 1.sam > possorted_bam_chrM_new.bam
samtools index possorted_bam_chrM_new.bam
rm 1.sam
samtools index possorted_bam_chrM.bam
java -Xmx8g -jar /public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar -R /public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/refdata-cellranger-atac-hg19-1.2.0/fasta/genome.fa -T RealignerTargetCreator -known /public/home/chenbzh5/DB/hg19/SNP/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /public/home/chenbzh5/DB/hg19/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -nt 8 -I possorted_bam_chrM.bam -o possorted_bam_chrM.target.intervals
java -Xmx10g -jar /public/home/chenbzh5/Tools/GenomeAnalysisTK_3.5-0.jar -R /public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/refdata-cellranger-atac-hg19-1.2.0/fasta/genome.fa -T IndelRealigner -maxReads 10000000 -known /public/home/chenbzh5/DB/hg19/SNP/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known /public/home/chenbzh5/DB/hg19/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf -I possorted_bam_chrM.bam -targetIntervals possorted_bam_chrM.target.intervals -o possorted_bam_chrM_realign.bam
samtools sort -t CR possorted_bam_chrM_realign.bam -o possorted_bam_chrM_realign_sorted.bam

