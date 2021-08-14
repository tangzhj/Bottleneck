human_coding_sequence: 
complete coding sequence in human mitochondral genome from ncbi
13 coding sequence in total
coding_clean——coding sequence without headers

human_Protein_sequence:
complete protein sequence encoded by human mitochondrial CDS
protein_clean——protein sequence without headers

human_same_protein:
file contains two same human protein sequence 
each sequence is a combination of 13 proteins encoded by human mitochondral CDS
as input file for clustalw protein alignment

human_same_protein.aln:
alignment result from clustalw

human_same_condon:
feasible input (sequence) file for PAML
change manually from human_coding_sequence.txt
because some CDS in human_coding_sequence.txt contain terminators or incomplete condons at the ends of sequences
also, human_same_condon.txt can be used as input file together with human_same_protein for PAL2NAL web edition
PAL2NAL link: http://www.bork.embl.de/pal2nal/index.cgi?example=Yes#RunP2N

human-mouse: compare ND1 gene in human and mouse mitochondrial genome
human/mouse_sequence.fasta: complete mitochondrial genome in human/mouse
two_sequence: a file contains both human and mouse complete mt genome
two_ND1_protein: a file contains both human and mouse ND1 gene sequence
two.tree: tree file——necessary input for PAML
human_mouse.ctl: control file——necessary input for PAML, determines input, output, model and parameters we use

Basic protocol to estimate N, S and dN/dS from fasta file:
(1)protein alignment: 
clustalw2 two_ND1_protein.txt

(2)convert original sequence to CDS: 
pal2nal.pl two_ND1_protein.aln two_sequence.txt -codontable 2 -output paml -nogap  >  two_condon.txt

(3)calculate N, S, dN/dS:
codeml human_mouse.ctl

Things need to be paid attention for:
(1)Since we can't combine 13 protein sequences together and extract CDS, we change human_coding_sequence.txt into a feasible input file for PAML
(因为coding sequence本身是完整的，而不同的CDS之间却会有间隔，PAL2NAL需要通过连续的氨基酸序列比对回相应密码子，所以并不能同时读取分割开的多个CDS)
(2)N, S result PAML prints to terminal may be different from result in output, why?
(3)Usually, PAML is used to estimate evolution relations between different species, but if you insist taking two same sequences as input, N and S can still be calculated, but not dN/dS
(N和S只是根据物种的DNA和氨基酸序列得到的，DNA序列上所有可能产生的非同义突变和同义突变的数目）

additional_test_for_single_gene:
在final中我们是将线粒体基因组上13个CDS（删去终止子和末端多余碱基）拼在一起作为输入计算N和S，但是N/S比例似乎比预想要低
所以我们将13个CDS分开各自计算了一遍N，S，发现分开计算的N/S比和拼起来计算的N/S比之间存在差异，目前我们使用的是分开计算的比值