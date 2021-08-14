library(dplyr)
# 1.germline filter
snv_merge = read.table(file = "/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/mk_realign/pool_output/CD34+1_realign_MPP/CD34+1_realign_MPPcell.merge.snv",
                       header = TRUE)
germline <- filter(snv_merge, 
                   Reads2/(Reads1 + Reads2) > 0.9)
                     # (Reads1 + Reads2) > mean(Reads1 + Reads2)*0.5)
write.table(germline, file = "./MPP_cell_realign/MPP_cell_realign.merge.q20Q30.germline.snv",
            quote = FALSE,
            sep = "\t",
            row.names = F)
somatic <- filter(snv_merge, 
                  Reads2/(Reads1 + Reads2) < 0.9 &
                    Reads2/(Reads1 + Reads2) > 0.05 &
                    (Reads1 + Reads2) > mean(Reads1 + Reads2)*0.5)
write.table(somatic,
            file = "./MPP_cell_realign/MPP_cell_realign.merge.q20Q30.somatic.snv",
            sep = "\t",
            quote = FALSE,
            row.names = F)

print("germline_filter")

# 2.sc_SNV filter
# define file path and an empty dataframe
arg_1 = "/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/mk_realign/CD34_rep1/realign2_MPP/"
files <- dir(arg_1)
path <- arg_1
i <- 1
final <- data.frame("Chrom" = 0 , "Position" = 0, "Ref" = 0, "VarAllele" = 0)
final <- final[-length(final$Chrom),]
sc_germline <- final

for(i in 1 : length(files)) {
  x <- paste0(path, files[i], "/", files[i],".snv")
  y <- read.table(file = x, header = T, colClasses = c("character"))
  y$Reads1 <- as.numeric(y$Reads1)
  y$Reads2 <- as.numeric(y$Reads2)
  y$Reads2Plus <- as.numeric(y$Reads2Plus)
  y$Reads2Minus <- as.numeric(y$Reads2Minus)
  z <- filter(y, Reads2Plus > 1 & 
                Reads2Minus > 1 &
                Reads2Plus/(Reads2Minus + Reads2Plus) < 0.7 &
                Reads2Plus/(Reads2Minus + Reads2Plus) > 0.3 &
                (Reads1 + Reads2) >= 20 &
                Reads2/(Reads1 + Reads2) >= 0.10 &
                !(y$Ref == "G" & y$VarAllele == "T") &
                !(y$Ref == "C" & y$VarAllele == "A")
  )
  k <- filter(z, Reads2/(Reads1 + Reads2) >= 0.90)
  n <- select(z, c("Chrom", "Position", "Ref", "VarAllele"))
  k <- select(k, c("Chrom", "Position", "Ref", "VarAllele"))
  final <- merge(final, n, all = T)
  # final <- rbind(final, n)
  sc_germline <- rbind(sc_germline, k)
  # sc_germline <- merge(sc_germline, k, all = T)
}

frequency <- as.data.frame(table(sc_germline$Position))
frequency <- filter(frequency, frequency$Freq > length(files)*0.9)
final <- unique(final)
final_remove <- c()
for (i in 1:length(final$Position)) {
  if (final$Position[i] %in% frequency$Var1){
    final_remove <- append(final_remove, i)
  }
}
if(length(final_remove >= 1)){
  final <- final[-final_remove,]
}

write.table(final, file = "./MPP_cell_realign/MPP_cell_realign.sc.filter.snv", sep = "\t", quote = F, row.names = F)

# input germline and blacklist
SNV_filter <- read.table(file = "./MPP_cell_realign/MPP_cell_realign.sc.filter.snv", header = T)
germline <- read.table(file = "./MPP_cell_realign/MPP_cell_realign.merge.q20Q30.germline.snv", header = T)
blacklist <- read.table(file = "./blacklist.txt", header = T)

# remove germline mutation
SNV_remove <- c()
for(i in 1:length(SNV_filter$Position)) {
  if(SNV_filter$Position[i] %in% germline$Position) {
    SNV_remove <- append(SNV_remove, i)
  }
}
if(length(SNV_remove) >= 1){
  SNV_filter <- SNV_filter[-SNV_remove,]}

# remove mutation in black list
blacklist_remove <-c()
for(i in 1:length(SNV_filter$Position)) {
  if(SNV_filter$Position[i] %in% blacklist$Position) {
    blacklist_remove <- append(blacklist_remove, i)
  }
}
if(length(blacklist_remove) >= 1){
  SNV_filter <- SNV_filter[-blacklist_remove,]}

# arrange
SNV_filter <- arrange(SNV_filter, Position)

# output
write.table(SNV_filter, file = "./MPP_cell_realign/MPP_cell_realign.sc.filter.rmgermline.rmblack.snv", sep = "\t", quote = F, row.names = F)
warnings()
print("sc_SNV_filter")

# 3.make SNV matrix
arg_1 <- "/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/mk_realign/CD34_rep1/realign2_MPP/"
arg_2 <- "./MPP_cell_realign/MPP_cell_realign.sc.filter.rmgermline.rmblack.snv"

# input and path
files = dir(arg_1)
path = arg_1
SNV_conf <- read.table(file = arg_2, header = T)

# create an empty matrix
SNV_matrix <- matrix(nrow = length(files),ncol = (length(SNV_conf$Position) + 2))
fre_matrix <- matrix(nrow = length(files), ncol = (length(SNV_conf$Position) + 2))

# fill SNV number and frequency matrix
for (t in 1:length(files)) {
  sc_SNV <- read.table(file = paste0(path, files[t], "/", files[t], ".counts"), colClasses = "character", header = F)
  sc_SNV$V9 <- as.numeric(sc_SNV$V9)
  sc_SNV$V4 <- as.numeric(sc_SNV$V4)
  SNV_matrix[t, 1] <- files[t]
  fre_matrix[t, 1] <- files[t]
  SNV_matrix[t, length(SNV_conf$Position) + 2] <- (sum(sc_SNV$V4)/length(sc_SNV$V4))
  fre_matrix[t, length(SNV_conf$Position) + 2] <- (sum(sc_SNV$V4)/length(sc_SNV$V4))
  for (i in 1:length(SNV_conf$Position)) {
    if (SNV_conf$Position[i] %in% sc_SNV$V2) {
      k <- which(sc_SNV$V2 == SNV_conf$Position[i])
      if (sc_SNV$V3[k] == SNV_conf$Ref[i] &
          sc_SNV$V6[k] == SNV_conf$VarAllele[i] &
          sc_SNV$V4[k] >= 8){
        SNV_matrix[t, i + 1] = sc_SNV$V9[k]
        fre_matrix[t, i + 1] = (sc_SNV$V9[k]/sc_SNV$V4[k])}
      else {SNV_matrix[t, i + 1] = 0
      fre_matrix[t, i + 1] = 0}
    }
    else {SNV_matrix[t, i + 1] = "na"
    fre_matrix[t, i + 1] = "na"}}
}

warnings()  
print("make_SNV_matrix")

# organize site information
site_infor <- c()
for (p in 1:length(SNV_conf$Position)){
  site_infor <- append(site_infor, paste0("chrM: ", SNV_conf$Position[p], ": ", SNV_conf$Ref[p], "/", SNV_conf$VarAllele[p]))
}

# output
write.table(SNV_matrix, file = "./MPP_cell_realign/SNV_matrix.txt", sep = "\t", quote = F, row.names = F)
write.table(fre_matrix, file = "./MPP_cell_realign/fre_matrix.txt", sep = "\t", quote = F, row.names = F)
write.table(site_infor, file = "./MPP_cell_realign/site_infor.txt", col.names = F, row.names = F, quote = F)

# 4.test for contamination
my_SNV <- read.table(file = "./MPP_cell_realign/site_infor.txt", header = F, sep = ":")
reference <- read.table(file = "/public/home/chenbzh5/DB/hg19/SNP/dbsnp_138.hg19.chrM.vcf", header = F)

# make empty dataframe
SNV_marked <- my_SNV[1,]
SNV_marked <- SNV_marked[-1,]

# match
for (i in 1:length(my_SNV$V2)){
  if (my_SNV$V2[i] %in% reference$V2){
    SNV_marked <- merge(SNV_marked, my_SNV[i,], all = T)
  }
  else{i <- i + 1}
}

write.table(SNV_marked, file = "./MPP_cell_realign/marked_site.txt", sep = "\t", quote = F, row.names = F, col.names = F)
# write.table(SNV_marked, file = "C:/Users/E/Desktop/discover mito DNA bottleneck in T cell/data and code/Follow_up_analysis/MPP_cell_realign/SNV_marked.txt", sep = "\t", colnames = F)

# overlapped mutation in single cell
path <- "/public/home/chenbzh5/project/mitoDNA_bottleneck/discover_bottleneck/output_GSE142745/PBMC_scATAC_1/mk_realign/CD34_rep1/realign2_MPP/"
files <- dir(path)
reference <- read.table(file = "/public/home/chenbzh5/DB/hg19/SNP/dbsnp_138.hg19.chrM.vcf", header = F)

# make an empty vector
cell <- files
contamination_rate <- c()

for (i in 1:length(files)){
  full <- 0
  marked <- 0
  filename <- paste0(path, files[i], "/", files[i], ".counts")
  sc_SNV <- read.table(file = filename, header = F)
  for (t in 1:length(sc_SNV$V2)) {
    if (sc_SNV$V2[t] %in% my_SNV$V2){
      full <- full + 1
      if (sc_SNV$V2[t] %in% reference$V2){
        marked <- marked + 1 
      }
    }
    else{t <- t + 1}
  }
  contamination_rate[i] <- round(marked/full, digits = 3)
}

total <- data.frame(cell, contamination_rate)
write.table(total, file = "./MPP_cell_realign/sc_contamination.txt", sep = "\t", quote = F, row.names = F, col.names = F)
