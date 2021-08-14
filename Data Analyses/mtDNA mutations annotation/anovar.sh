for file in `ls *.filter`; 
do awk '{print "chrM\t"$2"\t"$3"\t"$4"\t"$5}'  $file  > $file.format; 
/public/home/jinxu/software/annovar/table_annovar.pl $file.format  /public/home/jinxu/software/annovar/humandb/ -buildver hg19 -out $file -protocol MT_ensGene -operation g -remove ; 
done
##Input file 5.Data for plotting/matrix/*Site_infor.txt. Result in 5.Data for plotting/annovar_result.