#PBS -N cellranger
#PBS -j oe
#PBS -q batch
#PBS -S /bin/sh
#PBS -l nodes=1:ppn=24
#PBS -l mem=30000m


cd /public/home/tangzj/PB_ATAC/GSE142745/cellrange_output/
sample="SRR11539033"
fastq="/public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/fatsq"
ref="/public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/refdata-cellranger-atac-hg19-1.2.0"
#source /public/home/tangzj/PB_ATAC/GSE142745/SRR10804559/refdata-cellranger-atac-hg19-1.2.0/sourceme.basih
/public/home/jinxu/software/cellranger-atac-1.1.0/cellranger-atac count	 --id=SRR11539033 \
	--fastqs=$fastq \
	--localcores=12 \
	--reference=$ref \
	--localmem=10 \
