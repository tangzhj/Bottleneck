for file in `ls *.gz`
do
~/software/FastQC/fastqc $file  
done
