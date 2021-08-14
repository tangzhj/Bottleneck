export PATH=$PATH:/public/home/chenbzh5/project/mitoDNA_bottleneck/src
input="/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/data/GSE96769/LMPP"
output="/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE96769/LMPP"
project_name="LMPP"
# make sample list from raw fastq folder 
# input paramters : fastq_folder	 project_name	output_file 
0_Make_SampleTable1.pl   $input  $project_name  sample.csv
# THe core pipeline which align each single cell and split into genome and mito genome.
# input parameters: sample sheet from the first step	output_folder 
1_scATAC_mito.sh  sample.csv  $output
# SNP and peak calling with merged single cell data after alignment.
# input parameters	output_folder from last step and name for the 
2_merge_cells1.sh   $output  $project_name
# Generate summary  table
3_Get_SummaryTable.sh  sample.csv $output 





Command to run pipeline: sh runpipeline.sh

A pbs file template: GSE107816_memory.pbs
we use PBS system to submit and run our work, so we need to write a pbs file first

Command to submit pbs file: qsub -l nodes=1:ppn=20 GSE107816_memory.pbs

Contents and some modification in file runpipeline.sh we need to be aware of:
0_Make_SampleTable1.pl —— we do some minor adjustment
we use 0_Make_SampleTable1.pl in pipeline

1_scATAC_mito.sh —— with trimming process
1_scATAC_mito1.sh  —— without trimming process
choose the right file before we start pipeline

2_merge_cells1.sh —— set up the right path for software we need
so we use 2_merge_cells1.sh in pipeline instead of 2_merge_cells.sh

3_Get_SummaryTable.sh —— commands to make summary.tsv file
3_Get_SummaryTable.sh works well in pipeline
but if we run merge_cell twice or more for some reason, we'd better use 3_Get_SummaryTable1.sh

Specific commands for mapping and snv calling are recorded in file scATAC_mtSMC_commands.sh