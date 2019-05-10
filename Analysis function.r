# Mapping differential methylation sites to gene


###################################################################
# Find each cancer DMP corresponding to the promoter region, 
# and write to the file, and then from the genetic level
##################################################################

fileData <- list(
  "/home/songy/code_dmr/data_BLCA",
  "/home/songy/code_dmr/data_BRCA",
  "/home/songy/code_dmr/data_COAD",
  "/home/songy/code_dmr/data_ESCA",
  "/home/songy/code_dmr/data_HNSC",
  "/home/songy/code_dmr/data_KIRC",
  "/home/songy/code_dmr/data_KIRP",
  "/home/songy/code_dmr/data_LIHC",
  "/home/songy/code_dmr/data_LUAD",
  "/home/songy/code_dmr/data_LUSC",
  "/home/songy/code_dmr/data_PRAD",
  "/home/songy/code_dmr/data_THCA",
  "/home/songy/code_dmr/data_UCEC")

fileprobe <- "/home/songy/code_dmr/geneTSS/data/probe_full_annotation.txt"
probeData <- read.table(fileprobe, stringsAsFactors = F)  # locus info #ID  CHR	POS	QUALITY:GENE:CGI:SPEC:CpGs:DESIGN:TYPE:RANDOM
rownames(probeData) <- probeData[,1] 
# gene info
standardGeneData <- read.csv("/home/songy/code_dmr/geneTSS/data/standardGeneData.csv", stringsAsFactors = F)

file_gene <- "/home/songy/code_dmr/geneTSS"  # gene info file
sec_file <- "final_DMR_notCpG_mcfs"  # CpG island or non CpG island
#sec_file <- "final_DMR_notCpG_mcfs"
#onefile <- "/home/songy/code_dmr/data_COAD"
#x <- FromDMPGetGene(onefile,sec_file,file_gene,probeData,standardGeneData)

allfile <- lapply(fileData, FromDMPGetGene,sec_file,file_gene,probeData,standardGeneData)




################################################################
# The commonalities between cancers were found at the DMP level, 
# and then mapped to genes to see the characteristics of pathway
################################################################
fileData <- list(
  "/home/songy/code_dmr/data_BLCA",
  "/home/songy/code_dmr/data_BRCA",
  "/home/songy/code_dmr/data_COAD",
  "/home/songy/code_dmr/data_ESCA",
  "/home/songy/code_dmr/data_HNSC",
  "/home/songy/code_dmr/data_KIRC",
  "/home/songy/code_dmr/data_KIRP",
  "/home/songy/code_dmr/data_LIHC",
  "/home/songy/code_dmr/data_LUAD",
  "/home/songy/code_dmr/data_LUSC",
  "/home/songy/code_dmr/data_PRAD",
  "/home/songy/code_dmr/data_THCA",
  "/home/songy/code_dmr/data_UCEC")

fileprobe <- "/home/songy/code_dmr/geneTSS/data/probe_full_annotation.txt"
probeData <- read.table(fileprobe, stringsAsFactors = F)  # loci info
rownames(probeData) <- probeData[,1] 

standardGeneData <- read.csv("/home/songy/code_dmr/geneTSS/data/standardGeneData.csv", stringsAsFactors = F)

file_gene <- "/home/songy/code_dmr/geneTSS"  
sec_file <- "final_DMR_CpG_mcfs"  
times <- 10  

dmpInf <- GetCombinDMR(times, fileData, sec_file)
fileout_gene <- paste(file_gene, paste0("geneTSS_",dmpInf$tss,"_com_freq_",times,".txt"),sep = "/")
fileout <- DMPtoGene(dmpInf$site_Confirmed, probeData, standardGeneData, fileout_gene)


########################Differentially methylated loci were found in multiple cancers########################
allTimes <- 1:13
sec_file <- "ensemble"
a <- lapply(allTimes, GetCombinDMR,fileData, sec_file,cutoff)


fileData <- list(
  "/home/songy/code_dmr/data_BLCA",
  "/home/songy/code_dmr/data_BRCA",
  "/home/songy/code_dmr/data_COAD",
  "/home/songy/code_dmr/data_ESCA",
  "/home/songy/code_dmr/data_HNSC",
  "/home/songy/code_dmr/data_KIRC",
  "/home/songy/code_dmr/data_KIRP",
  "/home/songy/code_dmr/data_LIHC",
  "/home/songy/code_dmr/data_LUAD",
  "/home/songy/code_dmr/data_LUSC",
  "/home/songy/code_dmr/data_PRAD",
  "/home/songy/code_dmr/data_THCA",
  "/home/songy/code_dmr/data_UCEC")

GetSampleNumber <- function(filedata){
  load(file.path(filedata,"trainData.Rdata"))
  num <- dim(trainData)[2]
  
  cat(filedata, " ",num,"\n" )
  return(num)
}

num <- lapply(fileData,GetSampleNumber)





# fileData <- "F:/ChAMP_data_DMRs/1new_code_LR/data/refFlat_hg19.txt"
# 
# data <- read.table(fileData)
# data <- data[,1:6]
# 
# StandardGeneData <- function(x){
# #   str(x)
# #   cat(x)
#   if(x[4] == "+"){
#     a <- as.numeric(x[5])
#     geneD = c(a-2500, a+500)
#   }
#   else{
#     a <- as.numeric(x[6])
#     geneD = c(a-500, a + 2500)
#   }
#   # cat("sadf:",geneD, "\n")
#   return(c(x[1], x[2], x[3], geneD))
#     
# }
# standardGeneData <- apply(data, 1, StandardGeneData)
# standardGeneData <- t(standardGeneData)
# colnames(standardGeneData) <- c("gene_name", "gene_info", "chr", "TSS_start", "TSS_end")
# standardGeneData <- as.data.frame(standardGeneData)
# # standardGeneData$TSS_start <- as.numeric(standardGeneData$TSS_start)
# # standardGeneData$TSS_end <- as.numeric(standardGeneData$TSS_end)
# 
# fileW <- "F:/ChAMP_data_DMRs/1new_code_LR/data/standardGeneData.csv"
# write.csv(standardGeneData, file = fileW, row.names = F, quote = F)


SureDMRtoGene <- function(oneDMR, standardGeneData, probeData){

  
  probe <- probeData[oneDMR,]
  probe <- as.matrix(probe)
  site_num <- as.numeric(probe[1,3])
  gene <- subset(standardGeneData, chr == as.character(probe[1,2]) & TSS_start <= site_num & TSS_end >= site_num)
  
  if(dim(gene)[1] < 1){
    
    gene_name <- c("-1")
    gene_info <- c("-1")
    chr <- c("-1")
    TSS_start <- c(-1)
    TSS_end <- c(-1)
    site <- rep(as.character(probe[1,1]), 1)
    site_num <- rep(site_num, 1)
    gene <- data.frame(gene_name, gene_info, chr, TSS_start, TSS_end, site, site_num)
  }
  else{
    gene$site <- rep(as.character(probe[1,1]), dim(gene)[1])
    gene$site_num <- rep(site_num, dim(gene)[1])
  }
  
  return(gene)
}

DMPtoGene <- function(comDMR, probeData, standardGeneData,fileout_gene){

  
  gene_site <- lapply(comDMR, SureDMRtoGene, standardGeneData, probeData)
  gene_site <- Reduce(rbind, gene_site)
  name <- c("site", "site_num", "gene_name", "gene_info", "chr", "TSS_start", "TSS_end")
  gene_site <- gene_site[name]
  
  xx <- gene_site[,3]
  k1 <- table(xx)
  a <- names(k1)
  a <- a[which(a!="-1")]
  cat("tss_number:", length(a), " ")
  a <- matrix(a, nrow = length(a), ncol = 1)
  a <- as.data.frame(a)
  #str(a)
  xxx <- xx[xx=="-1"]
  cat("not_tssNumber:", length(xxx),"\n")
  
  write.table(a, file = fileout_gene, row.names = F, col.names = F, quote = F,sep = "\t")
  
  return(fileout_gene)
}

FromDMPGetGene <- function(onefile,sec_file,file_gene,probeData,standardGeneData){
  
  
  dmpInf <- GetCanserDMP(onefile,sec_file)
  #str(dmpInf)
  
 
  canser <- strsplit(onefile,split = "_")[[1]]
  canser <- canser[length(canser)]
  cat(canser," ")
  
  fileout_gene <- paste(file_gene, paste0("geneTSS_",dmpInf$tss,"_",canser,".txt"),sep = "/")
  fileout <- DMPtoGene(dmpInf$site_Confirmed, probeData, standardGeneData,fileout_gene)
  
  return(fileout)
}


############## DMP Common ###########################
GetCanserDMP <- function(onefile,sec_file){

  
  namef <- load(paste(onefile, paste0(sec_file,".Rdata"), sep="/"))

  if(namef == "final_DMR_CpG_mcfs"){
    site_Confirmed <- final_DMR_CpG_mcfs
    tss <- "CpG"
  }
  if(namef == "final_DMR_notCpG_mcfs"){
    site_Confirmed <- final_DMR_notCpG_mcfs
    tss <- "notCpG"
  }
  

  if(namef == "final_DMR_gene_mcfs"){
    site_Confirmed <- final_DMR_gene_mcfs
    tss <- "gene"
  }
  if(namef == "final_DMR_notGene_mcfs"){
    site_Confirmed <- final_DMR_notGene_mcfs
    tss <- "notGene"
  }
  if(namef == "final_DMR_TSS_mcfs"){
    site_Confirmed <- final_DMR_TSS_mcfs
    tss <- "TSS"
  }
  
  return(list(site_Confirmed=site_Confirmed,tss=tss))
}


GetCombinDMR <- function(times, fileData, whichSec){
 
  
  all_fea <- lapply(fileData, GetCanserDMP, whichSec)
  area <- all_fea[[1]]$tss
  all_fea <- unlist(lapply(all_fea, function(x){x$site_Confirmed}))
  fea_table <- sort(table(all_fea),decreasing=T)
  select_fea <- fea_table[which(fea_table>=times)]
  select_fea <- names(select_fea)
  cat(area, " ",times, " ", length(select_fea),"\n")
  
  site_Confirmed=select_fea
  tss=area
  return(list(site_Confirmed=site_Confirmed,tss=tss))
}