##################################################################################
##ABC prepare parameter data frame
##################################################################################
library(data.table)
library(stringr)
library(abc)
library(ggplot2)
library(ggpubr)
dirs <- c("7861e-7","27631e-7","10641e-7","7861e-8","27631e-8","10641e-8")
for (mydir in dirs){
  message("Calculating ",mydir,"...")
  setwd(paste("~/Simulation/",mydir,"/",sep = ""))
  files <- list.files()
  d <- 0
  for (my.file in files){
    x <- read.delim(my.file,header =F)
    setDT(x)
    e1 <- data.table()
    e2 <- data.table()
    for(j in 1:nrow(x)){
      temp = x[j,]
      df1 <- data.table(
        numCell = as.numeric(temp$V1),
        alpha = as.numeric(temp$V2),
        mu = as.numeric(temp$V3),
        Td = as.numeric(temp$V4),
        t = as.numeric(temp$V5),
        numMT = as.numeric(temp$V6))
      freqs = temp$V7
      freqs = gsub("\\[", "", freqs)
      freqs = gsub("\\]", "", freqs)
      freqs =as.numeric(strsplit(freqs, ",")[[1]])
      br = seq(0,1,by=0.05)
      ranges = paste(head(br,-1), br[-1], sep="-")
      freq   = hist(freqs, breaks=br, right=TRUE, plot=FALSE)
      freq = freq$counts#/70
      freq <- setDT(as.list(freq))[]
      setnames(freq, colnames(freq), ranges)
      e1 <- rbind(e1, df1)
      e2 <- rbind(e2, freq)
    }
    e <- cbind(e1,e2)
    d <- d+nrow(e)
    print(d)
    write.table(e, file = paste("~/Dilution/ABC/",mydir,"/ABC_",my.file, sep = ""), sep = "\t", row.names = F, quote = F)
  }
}
##merge abc files
for (mydir in dirs){
  message("Calculating ",mydir,"...")
  setwd(paste("E:/Dilution/ABC/",mydir,"/",sep = ""))
  files <- list.files()
  d <- data.table()
  for (my.file in files){
    x <- read.delim(my.file,header =T)
    x$numMT = round(500*x$alpha^x$Td)
    d <- rbind(d, x)
  }
  print(nrow(d))
  write.table(d[1:1000000,], file = paste("~/Dilution/ABC_",mydir,".txt", sep = ""), sep = "\t", row.names = F, quote = F)
}
##################################################################################
##abc parameter inference
##################################################################################
samples <- list(c("7861e-7","B"),c("7861e-8","B"),c("10641e-7","NK"),c("10641e-8","NK"),c("27631e-7","T"),c("27631e-8","T"))
pdf(file = "ABC_posterior.pdf", width = 8, height = 6)
res <- data.table()
for(my.sample in samples){
  message("Calculating ",my.sample[1],"_",my.sample[2],"...")
  d <- read.table(paste("ABC_",my.sample[1],".txt",sep = ""), header = T)
  d1 <- d[,c(2,4,5,6)]
  d2 <- d[, 8:ncol(d)]
  setDT(d1)
  setDT(d2)
  targ <- read.delim(paste("~/DataForSimul0418/",my.sample[2],"_mutation_freqencies_0418.txt",sep = ""),header =T)
  setDT(targ)
  targ <- targ[freq > 0.05]
  br = seq(0.05,1,by=0.05)
  ranges = paste(head(br,-1), br[-1], sep="-")
  freq   = hist(targ$freq, breaks=br, right=TRUE, plot=FALSE)
  freq = freq$counts#/sampleSize[my.sample]#*1647
  freq <- setDT(as.list(freq))[]
  setnames(freq, colnames(freq), ranges)
  ## ---------------------------------------------------------------------------
  ## method = rejection
  rej1 <- abc(target=freq, param=d1[,1:3], sumstat=d2, tol=.001, 
              method ="rejection")
  sum1 <- summary(rej1)
  x1 <- data.frame(rej1$unadj.values)
  x1$numMT = round(500*x1$alpha^(x1$Td-1))
  ## 
  sum1 <- as.data.frame.matrix(sum1)
  sum1$numMT <- round(500*sum1$alpha^(sum1$Td-1))
  sum1$var <- row.names(sum1)
  sum1$Sample = paste0(my.sample[1],"_",my.sample[2])
  sum1$Method = "rejection"
  res <- rbind(res,sum1)
  p1 <- ggplot(x1)+geom_density(aes(x=alpha))+geom_rug(aes(x=alpha,y=0),position = position_jitter(height = 0))+
    xlim(0,1)+
    theme_classic()+ggtitle(paste("Rejection_Posterior distribution of alpha",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggplot(x1)+geom_density(aes(x=Td))+geom_rug(aes(x=Td,y=0),position = position_jitter(height = 0))+
    xlim(0,30)+
    theme_classic()+ggtitle(paste("Posterior distribution of Td",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  p3 <- ggplot(x1)+geom_density(aes(x=t))+geom_rug(aes(x=t,y=0),position = position_jitter(height = 0))+
    xlim(10,40)+
    theme_classic()+ggtitle(paste("Posterior distribution of t",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  p4 <- ggplot(x1)+geom_density(aes(x=numMT))+geom_rug(aes(x=numMT,y=0),position = position_jitter(height = 0))+
    xlim(0,500)+
    theme_classic()+ggtitle(paste("Posterior distribution of numMT",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  print(ggarrange(p1,p2,p3,p4,ncol=2,nrow=2))
  ## ---------------------------------------------------------------------------
  ## method = neuralnet
  rej2 <- abc(target=freq, param=d1[,1:3], sumstat=d2, tol=.001, 
              method ="neuralnet", transf = c("none","none","none"),
              kernal ="epanechnikov")
  sum2 <- summary(rej2)
  x2 <- data.frame(rej2$adj.values)
  x2$numMT = round(500*x2$alpha^(x2$Td-1))
  sum2 <- as.data.frame.matrix(sum2)
  sum2$numMT <- round(500*sum2$alpha^(sum2$Td-1))
  sum2$Td <- round(sum2$Td)
  sum2$t <- round(sum2$t)
  sum2$var <- row.names(sum2)
  sum2$Sample = paste0(my.sample[1],"_",my.sample[2])
  sum2$Method = "neuralnet"
  res <- rbind(res,sum2)
  p1 <- ggplot(x2)+geom_density(aes(x=alpha))+geom_rug(aes(x=alpha,y=0),position = position_jitter(height = 0))+
    xlim(0,1)+
    theme_classic()+ggtitle(paste("Neuralnet_Posterior distribution of alpha",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggplot(x2)+geom_density(aes(x=Td))+geom_rug(aes(x=Td,y=0),position = position_jitter(height = 0))+
    xlim(0,30)+
    theme_classic()+ggtitle(paste("Posterior distribution of Td",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  p3 <- ggplot(x2)+geom_density(aes(x=t))+geom_rug(aes(x=t,y=0),position = position_jitter(height = 0))+
    xlim(10,40)+
    theme_classic()+ggtitle(paste("Posterior distribution of t",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  p4 <- ggplot(x2)+geom_density(aes(x=numMT))+geom_rug(aes(x=numMT,y=0),position = position_jitter(height = 0))+
    xlim(0,500)+
    theme_classic()+ggtitle(paste("Posterior distribution of numMT",my.sample[1],my.sample[2],sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))
  print(ggarrange(p1,p2,p3,p4,ncol=2,nrow=2))
}
dev.off()
res <- res[var %in% c("Weighted Median:","Median:")]
res[,numCell := c(rep(786,4),rep(1064,4),rep(2763,4))]
res[,N1 := c(rep(500,4),rep(300,4),rep(500,4))]
res[,mu := rep(c(1e-7,1e-8,1e-7,1e-8,1e-7,1e-8),each =2)]
write.table(res, file = "ABC_stats1.txt",sep = "\t", row.names = F,
            quote = F)
##################################################################################
##abc model selection
##################################################################################
pdf(file = "model_selection_m2.pdf", width = 8, height = 6)
samples <- list(c("7861e-7","B"),c("7861e-8","B"),c("27631e-7","T"),c("27631e-8","T"),c("10641e-7","NK"),c("10641e-8","NK"))
# samples <- list(c("786","B"),c("7861e-7","B"),c("1064","NK"),c("10641e-7","NK"),c("2763","T"),c("27631e-7","T"))
for(my.sample in 1:3){
  message("Calculating ",samples[my.sample*2-1][[1]][1]," ",samples[my.sample*2][[1]],"...")
  ## model 1
  d1 <- read.table(paste("E:/Dilution/ABC_",samples[my.sample*2-1][[1]][1],".txt",sep = ""), header = T)
  ## model 2
  d2 <- read.table(paste("E:/Dilution/ABC_",samples[my.sample*2][[1]][1],".txt",sep = ""), header = T)
  ##
  d <- rbind(d1, d2)
  ## model 1
  d1 <- d[,1:6]
  d2 <- d[, 8:ncol(d)]
  setDT(d1)
  setDT(d2)
  # message(targetSamples[my.sample])
  targ <- read.delim(paste("E:/OneDrive/2.mtDNA/DataForSimul0418/",samples[my.sample*2][[1]][2],"_mutation_freqencies_0418.txt",sep = ""),header =T)
  targ <- targ[targ$freq > 0.05,]
  br = seq(0.05,1,by=0.05)
  ranges = paste(head(br,-1), br[-1], sep="-")
  freq   = hist(targ$freq, breaks=br, right=TRUE, plot=FALSE)
  freq = freq$counts#/sampleSize[my.sample]#*1647
  freq <- setDT(as.list(freq))[]
  setnames(freq, colnames(freq), ranges)
  ## abc
  rej <- abc(target=freq, param=d1[,.SD, .SDcols = c("alpha","Td","t")], sumstat=d2, tol=.001, method ="rejection")
  s <- summary(rej)
  x <- data.frame(rej$unadj.values)
  #dist
  s <- which(rej$region == "TRUE")
  r1 <- rej$dist[s[s<=1000000]]
  r2 <- rej$dist[s[s>1000000]]
  dt <- rbind(data.table(dist = r1, mod = "1e-7"),
              data.table(dist = r2, mod = "1e-8"))
  print(table(dt$mod))
  print(wilcox.test(dt[mod == "1e-7",dist],dt[mod == "1e-8",dist])$p.value)
  print(ggplot(dt, aes(x=mod, y=dist, color = mod))+geom_boxplot()+theme_classic()+
    ggtitle(paste(samples[my.sample*2-1][[1]][2]," cell",sep = "_"))+
    theme(plot.title = element_text(hjust = 0.5))+stat_compare_means(label.x = 1.5))
}
dev.off()
##################################################################################
##plot examples
##################################################################################
inputFile = "simulation.example1.txt"
eg <- read.delim(inputFile,header = F)
setDT(eg)
params <- unique(eg$V1)
pdf(file = paste(substr(inputFile,12,19),"_1.pdf",sep = ""), width = 9, height = 12)
for(my.param in params){
  message("Calculating ", my.param,"...")
  x <- eg[eg$V1 == my.param,]
  # x <- x[101:200,]
  ## rejection
  df <- data.table()
  for(j in 1:nrow(x)){
    temp = x[j,]
    freqs = temp$V6
    freqs = gsub("\\[", "", freqs)
    freqs = gsub("\\]", "", freqs)
    freqs =as.numeric(strsplit(freqs, ",")[[1]])
    temp <- data.table(noSimul = as.character(j), freqs = freqs, mu = temp$V2)
    df <- rbind(df, temp)
  }
  df <- df[freqs >=0.01]
  ##targ
  targ <- read.delim(paste0(strsplit(my.param,"_")[[1]][2],"_mutation_freqencies_0418.txt"),header = T)
  setDT(targ)
  targ <- targ[freq >= 0.01]
  df$mu <- round(df$mu,2)
  d <- data.table()
  for(my.Simul in unique(df$noSimul)){
    temp <- df[noSimul == my.Simul]
    temp <- temp[sample(1:nrow(temp), nrow(targ), replace = TRUE)]
    d <- rbind(d, temp)
  }
  df <- d
  p1 <- ggplot(df[mu!=1], aes(x=freqs, color = noSimul))+geom_density()+
    labs(title=paste0(my.param,"_w/ bottleneck"),x="Frequency",y="Density")+theme_classic()+theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5))
  p2 <- ggplot(df[mu==1], aes(x=freqs, color = noSimul))+geom_density()+
    labs(title=paste0(my.param,"_w/o bottlenekc"),x="Frequency",y="Density")+theme_classic()+theme(legend.position="none")+
    theme(plot.title = element_text(hjust = 0.5))
  df$mu <- as.character(df$mu)
  # one simulation
  d <- data.table()
  for(my.mu in unique(df$mu)){
    temp <- df[mu==my.mu]
    temp <- temp[noSimul == sample(unique(temp$noSimul),1)]
    d <- rbind(d, temp)
  }
  m <- rbind(data.table(cat = d$mu,freq = d$freqs),data.table(cat = "empir",freq = targ$freq))
  p3 <- ggplot(df[mu==unique(df$mu)[1]], aes(x=freqs))+
    geom_histogram(position="identity",alpha = 0.5,bins=10)+
    labs(title=paste0(my.param,"_density"),x="Frequency")+theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))
  p4 <- ggplot(df[mu==unique(df$mu)[2]], aes(x=freqs))+
    geom_histogram(position="identity",alpha = 0.5,bins=10)+
    labs(title=paste0(my.param,"_density"),x="Frequency")+theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))
  
  print(ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2))
}
dev.off()
##################################################################################
##  mutation burdens in simulated data
##################################################################################
inputFile = "simulation.example1.txt"
eg <- read.delim(inputFile,header = F)
eg <- cbind(data.frame(do.call('rbind', strsplit(as.character(eg$V1),'_',fixed=TRUE))),eg[,2:ncol(eg)])
setDT(eg)
eg <- eg[X4 == "1e-07" & X3 == "rejection"]
eg$numCell <- gsub("1e-7","",eg$X1)
eg$numCell <- as.numeric(eg$numCell)
params <- unique(eg$X1)
pdf(file = paste(substr(inputFile,12,19),"_num_mutations.pdf",sep = ""), width = 9, height = 12)
df <- data.table()
for(my.param in params){
  message("Calculating ", my.param,"...")
  x <- eg[eg$X1 == my.param,]
  for(j in 1:nrow(x)){
    temp = x[j,]
    freqs = temp$V6
    freqs = gsub("\\[", "", freqs)
    freqs = gsub("\\]", "", freqs)
    freqs =as.numeric(strsplit(freqs, ",")[[1]])
    freqs <- freqs[freqs>=0.01]
    temp <- data.table(cat = temp$X2, noSimul = as.character(j), numMut = length(freqs)/temp$numCell, mu = temp$V2)
    temp[,mu:= ifelse(mu ==1, "W/O_dilution", "W_dilution")]
    df <- rbind(df, temp)
  }
  ##targ
}
p1 <- ggplot(df, aes(x=cat, y = numMut, color = mu))+geom_boxplot()+
  labs(title="Number of mutations",x="Lymphocytes",y="Mutation burden")+theme_classic()+#theme(legend.position="none")+
  theme(plot.title = element_text(hjust = 0.5))
print(ggarrange(p1,ncol = 2, nrow = 1))
dev.off()
# ##################################################################################
# ##
# ##################################################################################
