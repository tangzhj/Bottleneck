#"use*.recal.mpileup to get depth.count and count aver_depth >> stats.tsv"
spname="Treg_colon2" 
echo "$spname"
for file in `ls ./$spname\_result/ERR*/*.recal.mpileup`
do
array=(${file//// }) 
/public/home/tangzj/scRNG-seq/tools/pileup_inf_rj.pl $file > ./$spname\_result/${array[2]}/depth.count
aver_depth=`awk 'BEGIN{sum=0;num=1}''{sum+=$4;num++}''END{print sum/num}' ./$spname\_result/${array[2]}/depth.count`
echo "${array[2]}	$aver_depth " >> "${spname}_stats.tsv" 
done