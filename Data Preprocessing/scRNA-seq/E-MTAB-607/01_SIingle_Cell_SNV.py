import os
refenerce = "/public/home/tangzj/scRNG-seq/workflow/reference/"
input_path = "/public/home/tangzj/scRNG-seq/E-MTAB-6072/Tcm_colon2/"
state = 0
dirs = os.listdir(input_path)
lit = []
case = 3 #select step

for i in dirs:
	if "fastq.gz" in i: 
		lit.append(i.split("_")[0])
lit = set(lit)
print(lit,len(lit))

#step 1
if case == 1:#case==0 for STAR
    for i in lit:
	if os.path.exists("./Tcm_colon2_result/"+ i):
            print(i+" is exist")
	    continue
        os.system("mkdir ./Tcm_colon2_result/"+ i)
        fastq1 = input_path + i + "_1.fastq.gz"
        fastq2 = input_path + i + "_2.fastq.gz"
        spname = "./Tcm_colon2_result/"+i+"/"+i


        if os.system("STAR --genomeDir "+str(refenerce)+" --runThreadN 30 --readFilesIn "+fastq1+" "+fastq2+" --readFilesCommand zcat --outFileNamePrefix "+spname+" --outSAMtype BAM Unsorted  --outBAMsortingThreadN 20") == 0:
            print("STAR is done")
        else:
            print("STAR error")


        if os.system("java -jar /public/home/tangzj/scRNG-seq/tools/picard/picard.jar AddOrReplaceReadGroups I="+str(spname)+"Aligned.out.bam O="+str(spname)+"Aligned.out.AdGroup.bam RGID=4 RGLB=library1 RGPL=illumina RGPU=machine RGSM=sample") == 0:
            print("AddGroup is done")
        else:
            print("AddGroup error")




#step 2
if case == 2:#case==1 for racal.bam
    for i in lit:
	if os.path.exists("./Tcm_colon2_result/"+i+"/"+i+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.bam"):
	    print(i+" is exist")
	    continue
        #os.system("mkdir ./Tem_result/"+ i)
        fastq1 = input_path + i + "_1.fastq.gz"
        fastq2 = input_path + i + "_2.fastq.gz"
	spname = "./Tcm_colon2_result/"+i+"/"+i
        print("#####"+spname+"#####")
        if os.system("java -jar /public/home/tangzj/scRNG-seq/tools/picard/picard.jar SortSam I="+str(spname)+"Aligned.out.AdGroup.bam O="+str(spname)+"Aligned.out.AdGroup.sorted.bam SORT_ORDER=coordinate") == 0:
            print("SortSam is done")
        else:
            print("SortSam error")



        if os.system("java -jar /public/home/tangzj/scRNG-seq/tools/picard/picard.jar MarkDuplicates I="+str(spname)+"Aligned.out.AdGroup.sorted.bam O="+str(spname)+"Aligned.out.AdGroup.sorted.markdup.bam M=sampleAligned.out.AdGroup.sorted.metrics") == 0:
            print("MarkDuplicates is done")
        else:
            print("MarkDuplicates error")



        if os.system("gatk SplitNCigarReads -R /public/home/tangzj/scRNG-seq/workflow/human_g1k_v37/human_g1k_v37.fasta -I "+str(spname)+"Aligned.out.AdGroup.sorted.markdup.bam -O "+str(spname)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.bam") == 0:
            print("SplitNCigarReads is done")
        else:
            print("SplitNCigarReads error")




        if os.system("gatk BaseRecalibrator -R /public/home/tangzj/scRNG-seq/workflow/human_g1k_v37/human_g1k_v37.fasta -I "+str(spname)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.bam --use-original-qualities -O "+str(spname)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.bam.recal_table -known-sites /public/home/tangzj/scRNG-seq/workflow/reference/1000G_phase1.snps.high_confidence.hg19.sites.vcf -known-sites /public/home/tangzj/scRNG-seq/workflow/reference/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf") == 0:
            print("BaseRecalibrator is done")
        else:
            print("BaseRecalibrator error")



        if os.system("gatk ApplyBQSR -R /public/home/tangzj/scRNG-seq/workflow/human_g1k_v37/human_g1k_v37.fasta -I "+str(spname)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.bam -bqsr "+str(spname)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.bam.recal_table -O "+str(spname)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.bam") == 0:
            print("ApplyBQSR is done")
        else:
            print("ApplyBQSR error")

#step3
if case == 3:#case==3 for snv
    for i in lit:
        #os.system("mkdir ./Tem_result/"+ i)
        print("-------------------------"+i+"----------------------------")
        fastq1 = input_path + i + "_1.fastq.gz"
        fastq2 = input_path + i + "_2.fastq.gz"
        spname = "./Tcm_colon2_result/"+i+"/"+i
        if os.system("samtools mpileup -l /public/home/tangzj/scRNG-seq/workflow/reference/chrM.len -q 20 -Q 30 -f /public/home/tangzj/scRNG-seq/workflow/human_g1k_v37/human_g1k_v37.fasta  -x "+spname+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.bam > "+spname+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup") == 0:
            print("mpileup is done")
        else:
            print("mpileup error")
        
        if os.system("java -Xmx4g  -jar /public/home/tangzj/scRNG-seq/tools/VarScan.v2.3.9.jar  pileup2snp  "+spname+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup --min-var-freq 0.01  --min-reads2 2 > "+spname+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup.snv") == 0:
            print("snp is done")
        else:
            print("snp error")

