################################
A demo for preprocessing of data
################################
start <- Sys.time()

library(ChAMP)
library(stringr)

fileData <- "/home/songy/code_bishe/data_raw/data_LUAD"

source('/home/songy/code_bishe/code_DataPrep/DMRsData_preprocess.R')


fileInfo <- FindFileInfo(fileData)
fileSummyLung <- 'summary_lung.csv'  
fileDMRs <- paste(fileData, sep = "/","DNA_Methylation/JHU_USC__HumanMethylation450/Level_1")
filelLungTestSet <- SureEachIdatLabel_toIdataDir(fileInfo, fileSummyLung, fileData, fileDMRs)

myLoad=champ.load(directory = fileDMRs)
myNorm=champ.norm()
save(myLoad,file=file.path(fileData,"myLoad.Rdata"))
save(myNorm,file=file.path(fileData,"myNorm.Rdata"))

champ.SVD()
batchNorm=champ.runCombat()
limma=champ.MVP()
lasso=champ.lasso(fromFile =TRUE, limma=limma)



args <- commandArgs(TRUE)  
fileData <- args[1]  # "/home/scfan/Data/ForYing/code/data_KIRC"
# fileData <- "F:/ChAMP_data_DMRs/code_DMRs"  

start <- Sys.time()

library(ChAMP)
library(stringr)

#source('DMRsData_preprocess.R', encoding = 'UTF-8')
source('DMRsData_preprocess.R')


fileInfo <- FindFileInfo(fileData)
fileSummyLung <- 'summary_lung.csv'  
fileDMRs <- paste(fileData, sep = "/","DNA_Methylation/JHU_USC__HumanMethylation450/Level_1")
filelLungTestSet <- SureEachIdatLabel_toIdataDir(fileInfo, fileSummyLung, fileData, fileDMRs)

myLoad=champ.load(directory = fileDMRs)
myNorm=champ.norm()
save(myLoad,file=paste())

champ.SVD()
batchNorm=champ.runCombat()
limma=champ.MVP()
lasso=champ.lasso(fromFile =TRUE, limma=limma)



#filelLungTestSet <-  'F:/ChAMP_data_DMRs/code_DMRs/data_KIRC/KIRC_all_lung_test_set.csv'

namesOfBatch <- HowManyBatchOfData(filelLungTestSet, dropSmallBatch = 1)
standardBatch <- namesOfBatch[1:6]
otherBatch <- setdiff(namesOfBatch, standardBatch)

fileWriteInfo <- paste(fileData, sep = "/","DNA_Methylation/JHU_USC__HumanMethylation450/Level_1")
dealStandardBatchData <- DealSomeBatchData(NULL, standardBatch, filelLungTestSet,fileWriteInfo, fileData)
save(dealStandardBatchData, file=paste0(fileData,"/afterComBatData/dealStandardBatchData.Rdata"))
dealAllBatchData <- lapply(otherBatch, DealSomeBatchData,standardBatch,filelLungTestSet, fileWriteInfo, fileData)
save(dealAllBatchData, file=paste0(fileData,"/afterComBatData/dealAllBatchData.Rdata"))


# load("afterComBatData/dealStandardBatchData.Rdata")
# load("afterComBatData/dealAllBatchData.Rdata")
standardExample <- dealStandardBatchData$col.name  
com.siteOfAllData <- lapply(dealAllBatchData, function(x){return(x$row.name)})
com.site <- Reduce(intersect, com.siteOfAllData) 

#infoStandData <- ComSiteExample_stand(dealStandardBatchData, com.site)
infoStandData <- ComSiteExample_v1(dealStandardBatchData, standardExample, com.site, fileData)
save(infoStandData, file=paste0(fileData,"/finalDataOfBatch/infoStandData.Rdata"))
infoFinaData <- lapply(dealAllBatchData, ComSiteExample_v1, standardExample, com.site, fileData)
save(infoFinaData, file=paste0(fileData,"/finalDataOfBatch/infoFinaData.Rdata"))


# load("finalDataOfBatch/infoStandData.Rdata")
# load("finalDataOfBatch/infoFinaData.Rdata")

#base <- "finalDataOfBatch"
fileSplit <- strsplit(infoStandData[[1]], split="/")[[1]]
base <- fileSplit[length(fileSplit)-1]
fileBase <- paste(fileData, base, sep = "/")
# filePreproData <- UnionData(fileBase, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)
# 

# #filePreproData <- "Data_prepro.csv"
# fileTrainData <- ConvertPreproToModelData(filePreproData, fileData) 

fileTrainData <- UnionAndConvertDataToModel(fileBase, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE)


end <- Sys.time()
end - start

