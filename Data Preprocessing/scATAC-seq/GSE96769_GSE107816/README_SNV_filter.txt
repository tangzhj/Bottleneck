for_group_process——将输入输出路径写入input/output_path.txt，即可批量处理不同类型细胞的SNV筛选
for_single_cell_type——一次只处理一类细胞的SNV筛选

Basic filtering rules:

1. Remove germline mutation:
(1) mutation frequency > 0.9 in more than 90% cells
(2) mutation frequency > 0.9 in merge.snv file

2. Consider strand balance:
(1) the number of Reads2Plus and ReadsMinus should be similar 
(即同一位点从双端测到的reads数应该相近)

(2) 0.3 < Reads2Plus/(Reads2Plus + Reads2Minus)  <0.7

3. Remove certain kind of mutation:
(1) remove G ——> T
(2) remove C ——> A

4. Remove sites in blacklist

5. Restriction of minimum frequency and sequening depth:
(1) Frequency > 0.1
(2) depth (Reads1 + Reads2) > 20

Filtering process:
1. Identify germline mutation from file merge.snv (rules1-(1) frequency > 0.9)
从merge的SNV位点中选出germline mutation

2. Choose confident somatic mutations from file single_cell.snv (rules2, 3, 5)
根据链平衡，测序深度，突变频率以及突变类型从每个细胞中确定可信的SNV位点
 
3. remove germline mutation
根据第一步得到的位点，以及90%中都含有的频率大于0.9的位点，记为germline mutation，并从第二步得到的SNV中删除

4. remove blacklist 
删除blacklist中记录的不可信位点

5. make SNV matrix
构建可用于作图以及其他分析的SNV矩阵