# STEP 1: filter mtSNV
# $17>=1 && $18>=1 ( make sure both positive and negative strand were called 
# Minor  Allele frequency more than 0.01 
# remove G->T and C->A, "N",based on previous empirical evidence. 
# sequence depth >= 20 for each cell, the sequence depth criteria should be adjust according to the data. 
path="/public/home/tangzj/scRNG-seq/E-MTAB-6072/"
array=(Treg_colon2)
for name in ${array[@]}
do
Project_infor=$path$name
#delete formal file
`rm $Project_infor\_result/varition/*.filter`
for file in `ls  $Project_infor\_result/ERR*/*.snv`
do
awk '($17>1 && $18>1)' $file | sed '1,1d' | awk '$6/($6+$5)>=0.1' | awk '(0.3<=$17/($17+$18) && $17/($17+$18)<=0.7)' | awk '!($3=="G" && $19=="T")' | awk '!($3=="C" && $19=="A")' | awk '$3!="N"'  | awk ' $6+$5>=20' >> $name"_result/varition/$name.mtDNA.snv.filter"
done
#awk '{print $1"\t"$2"\t"$3"\t"$19}' $name"_result/varition/$name.mtDNA.snv.filter" | sort -k2n -u >$name"_result/varition/$name.mtDNA.snv.filter.cut.uniq"
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$17"\t"$18"\t"$19}' $name"_result/varition/$name.mtDNA.snv.filter" | sort -k2n -u >$name"_result/varition/$name.mtDNA.snv.filter.cut.uniq"
done ## if running new dataset this should be deleted


# following is creat merge.recal.mpileup.snv to get *snv.germline

#spname=$Project_infor
#samtools merge -f $spname\_result/varition/merge.recal.bam  $spname\_result/ERR*/*.recal.bam
#samtools index $spname\_result/varition/merge.recal.bam
#samtools mpileup  -l  /public/home/tangzj/scRNG-seq/workflow/reference/chrM.len -q 20  -Q 30 -f /public/home/tangzj/scRNG-seq/workflow/human_g1k_v37/human_g1k_v37.fasta -x $spname\_result/varition/merge.recal.bam > $spname\_result/varition/merge.recal.mpileup
#java -Xmx4g  -jar /public/home/tangzj/scRNG-seq/tools/VarScan.v2.3.9.jar  pileup2snp  $spname\_result/varition/merge.recal.mpileup --min-var-freq 0.001  --min-reads2 2 > $spname\_result/varition/merge.recal.mpileup.snv
#echo $Project_infor
#done  ##done should be replaced to new one for running new dataset.


#germline=$path."memory.merge.mt.mpileup.q20Q30.snv.germline"
#blacklist="backlist.txt"
# filter germline mutations and further remove mutations in blacklist 
#fetch.pl  -m  0 -a 3 -b 2  $germline $Project_infor.mtDNA.snv.filter.cut.uniq   > $Project_infor.mtDNA.snv.filter.cut.uniq.rmgermline
#fetch.pl  -m  0 -a 2 -b 2  $blacklist $Project_infor.mtDNA.snv.filter.cut.uniq.rmgermline  > $Project_infor.mtDNA.snv.filter.cut.uniq.rmgermline.rmblacklist

