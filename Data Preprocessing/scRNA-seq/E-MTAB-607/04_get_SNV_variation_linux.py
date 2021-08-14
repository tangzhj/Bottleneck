import os
#two dic. first floor is cell and merging position ,second is position and mt frequence
spname = ["Tem_colon2","Tcm_colon2","Treg_colon2"]
name = spname[2]
path = "/public/home/tangzj/scRNG-seq/E-MTAB-6072/"
Project_infor = path+name+"_result/"  #spname is a lit
def addtwodimdict(thedict, key_a, key_b, val):
  if key_a in thedict.keys():
    thedict[key_a].update({key_b: val})
#celllit is cellname,dic_cell is a dict.key is cellname value is all cell snv position
dirs = os.listdir(Project_infor)#dirs for snv result
dic_cell = {}
celllit = []
blit = []
longlit = []
for d in dirs:
    if "ERR" in d:
        celllit.append(d)
    if os.path.exists(Project_infor+str(d)+"/"+str(d)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup.snv"):
        with open(Project_infor+str(d)+"/"+str(d)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup.snv",'r') as f:
            lines = f.readlines()
            if len(lines) > 2:
                
                for i in lines:
                    i = i.strip()
                    if "VarAllele" not in i:
                        lit = i.split("\t")
                        position = (str(lit[0])+":"+str(lit[1])+":"+str(lit[2])+"/"+str(lit[18]))
                        longlit.append(position)#All cell snv position

for i in celllit:
    dic_cell[str(i)] = {}
    for j in longlit:
        dic_cell[i].update({str(j):0})                       
                
                
for d in dirs:
    if os.path.exists(Project_infor+str(d)+"/"+str(d)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup.snv"):
        with open(Project_infor+str(d)+"/"+str(d)+"Aligned.sortedByCoord.out.AdGroup.sorted.markdup.split.recal.mpileup.snv",'r') as f:
            lines = f.readlines()
            if len(lines)>2:
                for i in lines:
                    i = i.strip()
                    if "VarAllele" not in i:
                        lit = i.split("\t")
                        position = (str(lit[0])+":"+str(lit[1])+":"+str(lit[2])+"/"+str(lit[18]))
                        var = float(str(lit[6]).replace("%",""))/100.0
                        
                        addtwodimdict(dic_cell,d,position,var)
                            
                             
                             
ff = open("/public/home/tangzj/scRNG-seq/E-MTAB-6072/"+name+"_stats.tsv",'r')
a = ff.readlines()
dic_depth = {}
for i in a:
  if "SampleID" not in i:
    depth  = float(i.split("\t")[1])
    cell_name = str(i.split("\t")[0])
    if cell_name not in dic_depth.keys():
      dic_depth[str(cell_name)] = depth
ff.close()                             
                             
                             
#print(str(dic_cell['SRR6358327']['chrMT:11051:A/T']))                           
                             
if os.path.exists(Project_infor+"/varition/"+name+"_matrix_variation.txt"):
    os.system("rm "+Project_infor+"/varition/"+name+"_matrix_variation.txt")
else:
    pass                                                        
                            
f = open(Project_infor+"/varition/"+name+"_matrix_variation.txt","a")#this is add
for j in dic_cell.keys():
    title = " "
    for chrm in dic_cell[j].keys():
        title= title+"\t"+str(chrm)
    title = title + "\tdepth\n"
    f.write(title)
    break
for j in dic_cell.keys():
    one_cell = str(j)
    for k in dic_cell[j].keys():
        one_cell = one_cell+"\t"+str(dic_cell[j][k])
    one_cell = one_cell + "\t"+str(dic_depth[j])+"\n"
    f.write(one_cell)
f.close()

        
