#!/bin/python
import pysam
import os
import subprocess
# file to split on
unsplit_file = "atac_possorted_bam.sorted.bam"
# where to place output files
out_dir = "./splitbam/"

# variable to hold barcode index
CB_hold = 'unset'
itr = 0
# read in upsplit file and loop reads by line
samfile = pysam.AlignmentFile( unsplit_file, "rb")
lit = []
for read in samfile.fetch( until_eof=True):
    # barcode itr for current read
    try:
        CB_itr = read.get_tag('CB')
    except:
        continue
    CB_itr = read.get_tag('CB')
    # if change in barcode or first line; open new file  
    if( CB_itr!=CB_hold or itr==0):
        # close previous split file, only if not first read in file
        if( itr!=0):
            if len(lit)>1:
                split_file = pysam.AlignmentFile( out_dir + "{}.bam".format(CB_hold), "wb", template=samfile)
                for kk in lit:
                    split_file.write(kk)
                split_file.close()
                lit = []
            lit = []
        CB_hold = CB_itr
        itr+=1
       # print(CB_hold)

    lit.append(read)
   # print(len(lit))
        #split_file = pysam.AlignmentFile( out_dir + "CB_{}.bam".format(CB_hold), "wb", template=samfile)
    
split_file.close()

samfile.close()
