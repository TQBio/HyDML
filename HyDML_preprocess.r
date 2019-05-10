#----------------------------------------------------#
# Programming instructions：Data preprocessing - data conversion, standardization, batch effects
# ---------------------------------------------------#


FindFileInfo <- function(fileData){
  # Get the address of the file information---eg."jhu-usc.edu_LIHC.HumanMethylation450.1.13.0.sdrf.txt"
  # 
  # Agrs：
  # fileData：Data folder
  #
  # Resturn：
  # address of the file information
  
  SureFileInfo <- function(xtxt){
    split.x <- unlist(strsplit(xtxt,'[.]'))
    sdrf <- split.x[length(split.x)-1] 
    if(sdrf == "sdrf")
      return(xtxt)
    else
      return(NULL)
  }
  
  pattern = "txt$"
  ignore.case = TRUE
  recursive = TRUE
  verbose = FALSE
  base <- paste(fileData, "METADATA/JHU_USC__HumanMethylation450", sep = "/")
  if (!all(file.exists(base))) 
    stop("'base' does not exists")
  info <- file.info(base)
  if (!all(info$isdir) && !all(!info$isdir)) 
    stop("'base needs to be either directories or files")
  if (all(info$isdir)) {
    txtfiles <- list.files(base, recursive = recursive, pattern = pattern, 
                           ignore.case = ignore.case, full.names = TRUE)
    if (verbose) {
      cat("[read] Found the following txt files:\n")
      print(txtfiles)
    }
  }
  else txtfiles <- list.files(base, full.names = TRUE)
  
  fileTxt <- lapply(txtfiles, SureFileInfo)
  fileTxt <- unlist(fileTxt)
  return(fileTxt)
  
}


FindNotExistData <- function(dataInfo,fileDMRs){
  # Remove undetected loci according to the annotation file
  # 
  # Agrs：
  # dataInfo: Idat data extracted from the annotation file
  # fileDMRs: Data folder 
  #
  # Resturn：
  # data annotated in file
  
  PasteFiles <- function(x){
    info <- paste(x[2], x[3], "Grn.idat", sep="_")
    return(c(x[1], info))
  }
  
  filesGrn <- list.files(fileDMRs, recursive = T, pattern = "Grn.idat", 
                         ignore.case = T, full.names = F)
  
  
  dataInfo.new <- dataInfo[c('Sample_Name','Sentrix_ID','Sentrix_Position')]
  filesInfo <- apply(dataInfo.new,1, PasteFiles)
  notExistFiles <- setdiff(filesInfo[2,], filesGrn)
  cat("data not exist in DNA_Methylation:", notExistFiles, "\n")
  existFiles <- filesInfo[2,] %in% notExistFiles
  existExample <- filesInfo[,!existFiles][1,]
  existExampleData <- dataInfo[,'Sample_Name'] %in% existExample
  existDataInfo <- dataInfo[existExampleData,]
  
  return(existDataInfo)
}

SureEachIdatLabel <- function(fileInfo, fileSummyLung, fileData, fileDMRs){
  # KIRC, for example
  # Determine the category of each Idat raw data (cancer or normal)
  #           Batch information，writte into KIRC_all_lung_test_set.csv
  #
  # Args：
  # fileInfo <-  'METADATA/jhu-usc.edu_KIRC.HumanMethylation450.1.11.0.sdrf.txt'
  # fileSummyLung <- 'summary_lung.csv'  Standard data form
  # fileData：Data folder
  # fileDMRs：DNA_Methylation data folder
  #
  # Returns:
  #'KIRC_all_lung_test_set.csv'
  
  Info = readLines(fileInfo, n = -1)
  getData <- function(x){
    data = unlist(strsplit(x, split='\t'))
    return(data)
  }
  dataSet <- lapply(Info, getData)
  dataSet <- as.data.frame(dataSet)
  dataSet <- t(dataSet)
  Data <- dataSet[-1,]
  colnames(Data) <- dataSet[1,]
  rownames(Data) <- NULL
  
  # To extract information about whether a tissue is cancer or not
  useInfo <- Data[,c('Comment [TCGA Barcode]','Array Data File')]
  getSampleInfo <- function(x){
    label <- unlist(strsplit(x[1], split='-'))[4]
    label <- substr(label, 1, 2)
    if (label <= '09' & label >= '01'){
      Sample_Group <- 'T'
    }
    if (label <= '19' & label >= '10'){
      Sample_Group <- 'N'
    }
    if (label <= '29' & label >= '20'){
      Sample_Group <- 'C'
    }
    
    Sentrix_ID <- unlist(strsplit(x[2], split='_'))[1]
    Sentrix_Position <- unlist(strsplit(x[2], split='_'))[2]
    
    ret = c(Sample_Group, Sentrix_ID, Sentrix_Position)
  }
  
  newDataInfo <- apply(useInfo, 1, getSampleInfo)
  newDataInfo <- t(newDataInfo)
  newDataInfo <-  unique(newDataInfo) 
  colnames(newDataInfo) <- c('Sample_Group', 'Sentrix_ID', 'Sentrix_Position')
  
  TN_DataInfo <- newDataInfo[which(newDataInfo[,1]=='T' |newDataInfo[,1] =='N'),]
  TN_DataInfo <- as.data.frame(TN_DataInfo)
  
  #TN_DataInfo['Sample_Name'] <- TN_DataInfo['Sample_Group']
  name <- c()
  for(i in 1:nrow(TN_DataInfo)){
    i_name <- paste(TN_DataInfo[i, 'Sample_Group'], i, sep='')
    name <- c(name,i_name)
  }
  TN_DataInfo['Sample_Name'] <- name
  
  otherName <- c('Sample_Plate','Pool_ID','Project','Sample_Well')
  TN_DataInfo[otherName] <- NA

  allName <- c('Sample_Name','Sample_Plate','Sample_Group','Pool_ID','Project',
               'Sample_Well','Sentrix_ID','Sentrix_Position')
  finalDataInfo <- TN_DataInfo[allName]
  

  lung_test <- read.csv(fileSummyLung, header=F, stringsAsFactor=F)
  nessaryInfo <- lung_test[1:8,]
  colnames(nessaryInfo) <- c('Sample_Name','Sample_Plate','Sample_Group','Pool_ID','Project',
                             'Sample_Well','Sentrix_ID','Sentrix_Position')
  
  TN_DataInfo.exist <- FindNotExistData(TN_DataInfo,fileDMRs)
  lung_test_set <- rbind(nessaryInfo, TN_DataInfo.exist)
  write.table(lung_test_set,file=paste(fileData, 'KIRC_all_lung_test_set.csv',sep = "/"),na = "",
              row.names = FALSE, col.names = FALSE, quote = F,  sep = ",")
  
  return(paste(fileData, 'KIRC_all_lung_test_set.csv',sep = "/"))
  
}



HowManyBatchOfData <- function(fileBatch, dropSmallBatch){
  # Function description: find which batches are in the data and return the batch name
  #
  # Args：
  # fileBatch <- 'KIRC_all_lung_test_set.csv'
  # dropSmallBatch：Remove small batches
  #
  # Returns:
  # The batch name contained in the data
  
  data <- read.csv(fileBatch, skip = 7,stringsAsFactors=F)
  Sentrix_ID <- as.factor(data$Sentrix_ID)
  #message('The number of samples contained in different batches：',table(Sentrix_ID))
  sortTable <- sort(table(Sentrix_ID), decreasing = T)
  namesOfBatch <- names(sortTable)[which(sortTable > dropSmallBatch)]
  dropName <- setdiff(names(sortTable), namesOfBatch)
  cat("drop name of Batch\n")
  cat(dropName, "\n")
  namesOfBatch <- as.numeric(namesOfBatch)
  return(namesOfBatch)
}



GetOneBatchData <- function(x,data){
  # Function description: select data belonging to a batch
  #
  # Args:
  # x:batch number
  # data：
  #
  # Returns:
  # All batch data information
  return(subset(data,Sentrix_ID==x ))
}

PreproBatchData <- function(filelLungTestSet, batchName, fileWriteInfo){
  # KIRC, for example ：Writes data information for multiple batches to a file
  # PreproBatchData(filelLungTestSet, batchNamen, fileWriteInfo)
  #
  # Args:
  # filelLungTestSet <- 'KIRC_all_lung_test_set.csv' 
  # batchName <- c(6042324009,6055432001,6042324043) 
  # fileWriteInfo <- paste(getwd(),'DNA_Methylation')
  #
  # Returns:
  # The location of multiple batches of information
  summyData <- read.csv(filelLungTestSet, nrows = 8,header = F, stringsAsFactors=F)
  data <- read.csv(filelLungTestSet, skip = 7, stringsAsFactors=F)
  selectData <- lapply(batchName, GetOneBatchData, data)
  selectData.new <- do.call("rbind", selectData)
  colnames(selectData.new) <- colnames(summyData)
  lungDataSet <- rbind(summyData, selectData.new)
  write.table(lungDataSet,file=paste(fileWriteInfo, 'KIRC_lung_test_set.csv', sep = '/'), 
              na = "",row.names = FALSE, col.names = FALSE, quote = F,  sep = ",")
  
  return(fileWriteInfo)
  
}


PropreDataToKIRC <- function(testdir,one.batch,fileData){
  # data pre-processing: Standardization, Remove batch effect
  #
  # Args：
  # testdir = paste(getwd(),'DNA_Methylation', sep = "/")
  # one.batch：The name of a particular batch
  # fileData: Data location
  
  myLoad <- champ.load(directory = testdir)
  myNorm <- champ.norm(beta = myLoad$beta, rgSet = myLoad$rgSet,
                       pd = myLoad$pd, mset = myLoad$mset, norm = "SWAN")
  
  resultsDir <- paste(fileData, "afterComBatData",sep="/") 
  if(!file.exists(resultsDir))
  {
    dir.create(resultsDir)
    message("Creating results directory. Results will be saved in ", resultsDir)
  }
  
  if (is.na(min(myNorm$beta))){
    file.myLoad <- paste0(fileData, '/afterComBatData/myLoad_', one.batch, '.Rdata')
    file.myNorm <- paste0(fileData, 'afterComBatData/myNorm_', one.batch, '.Rdata')
    file.dropExamole <- peste0(fileData, '/afterComBatData/dropExamole.txt')
    save(myLoad, file = file.myLoad)
    save(myNorm, file = file.myNorm)
    notNA.myNorm <- na.omit(t(myNorm$beta))
    nameNA <- rownames(notNA.myNorm)
    nameNorm <- colnames(myNorm$beta)
    dropEaxmple.name <- setdiff(nameNorm, nameNA)
    dropEaxmple <- c(one.batch,dropEaxmple.name)
    write.table(matrix(dropEaxmple,nrow=1), file=file.dropExamole, append=TRUE, row.names=F,  col.names = F)
    
    pd <- myLoad$pd
    newName <- rownames(pd) %in% dropEaxmple.name
    pd.new <- pd[!newName,]  
    batchNorm <- champ.runCombat(beta.c = t(notNA.myNorm), pd = pd.new)
    return(batchNorm)

    # return(t(notNA.myNorm))
  }
  #champ.SVD()
  batchNorm <- champ.runCombat(beta.c = myNorm$beta, pd = myLoad$pd)
  return(batchNorm)
  
  # return(myNorm$beta)
}


DealSomeBatchData <- function(oneBatch, standardBatch, filelLungTestSet, fileWriteInfo, fileData){
  # Six batches were selected as standard batches, and one batch was added for data preprocessing
  # 
  # Args：
  # oneBatch：A batch
  # standardBatch：Standard batch
  # filelLungTestSet：'KIRC_all_lung_test_set.csv' ,file name
  # fileWriteInfo：Write the data information to the methylated data set
  # fileData:Data location
  #
  # Returns：
  # information of batchNorm
  
  someBatch <- append(standardBatch, oneBatch)
  testdir  <-PreproBatchData(filelLungTestSet, someBatch, fileWriteInfo)
  batchNorm <- PropreDataToKIRC(testdir,oneBatch, fileData)
  # batchNorm <- as.data.frame(batchNorm)  # 删除
  batchNorm <- as.data.frame(batchNorm$beta)
  fileBat <- paste0(fileData, '/afterComBatData/afterComBatData_',oneBatch,'.csv')
  write.csv(batchNorm,file=fileBat)
  return(list(info=c(fileBat,dim(batchNorm)), col.name=colnames(batchNorm),
              row.name=rownames(batchNorm)))
}

ComSiteExample_stand <- function(oneBatch, com.site, fileData){
  # Function description: for a CSV file after correcting batch effect, extract the same site data
  #
  # Args：
  # oneBatch：A batch CSV file information
  # com.site：Name of common site
  # fileData: Data location
  #
  # Returns:
  # list：fileInfo--csv file information，siteInfo--site information
  
  data <- read.csv(oneBatch$info[1], stringsAsFactors =F)
  rownames(data) <- data[,1]
  newData <- data[,-1]
  
  comNewData <- newData[com.site, ]
  comNewData <- t(comNewData)
  
  fileFialData <- paste0(fileData, "/DMRs_data.txt")
  write.table(comNewData, file = fileFialData, sep = ",")
  
  return(list(fileInfo=fileFialData, siteInfo=colnames(comNewData)))
}

ComSiteExample_v1 <- function(oneBatch, standardExample, com.site, fileData){
  # For a batch CSV file after correcting the batch effect, remove the standard batch and extract the same site data
  #
  # Args：
  # oneBatch：A batch CSV file information
  # standardExample：Standard batch name
  # com.site：Name of common site
  # fileData:Data location
  #
  # Returns:
  # list：fileInfo--csv file information，siteInfo--site information
  
  data <- read.csv(oneBatch$info[1], stringsAsFactors =F)
  rownames(data) <- data[,1]
  newData <- data[,-1]
  
  comNewData <- newData[com.site, ]
  if (dim(comNewData)[2] == length(standardExample)){
    dropStandData <- comNewData
  }
  else{
    colNewData <- colnames(comNewData) %in% standardExample
    dropStandData <- comNewData[!colNewData]
  }
  
  
  batchName <- str_extract(oneBatch$info[1], "\\d+")
  
  resultsDir <- paste(fileData, "finalDataOfBatch",sep="/")  
  if(!file.exists(resultsDir))
  {
    dir.create(resultsDir)
    message("Creating results directory. Results will be saved in ", resultsDir)
  }		
  
  fileFialData <- paste0(fileData, "/finalDataOfBatch/finalData_", batchName, ".csv")
  write.csv(dropStandData, file = fileFialData)
  
  return(list(fileInfo=fileFialData, siteInfo=rownames(dropStandData)))
}


ComSiteExample_v2 <- function(oneBatch, standardExample, com.site, fileData){
  # Ditto
  #           
  #
  # Args：
  # oneBatch：
  # standardExample：
  # com.site：
  # fileData:
  #
  # Returns:
  # 
  
  data <- read.csv(oneBatch$info[1], stringsAsFactors =F)
  rownames(data) <- data[,1]
  newData <- data[,-1]
  
  comNewData <- newData[com.site, ]
  colNewData <- colnames(comNewData) %in% standardExample
  dropStandData <- comNewData[!colNewData]
  dropStandData <- t(dropStandData)
  
  fileFialData <- paste0(fileData, "/DMRs_data.txt")
  write.table(dropStandData, file = fileFialData, append = T, col.names = F, sep = ",")
  
  return(list(fileInfo=fileFialData, siteInfo=colnames(dropStandData)))
}


#base <- "E:/Rproject/ChAMP_data_DMRs/data_KIRC/finalDataOfBatch"
UnionData <- function(base, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE){
  # Combine data，rows represent samples ，columns represents loci
  #
  # Args：
  # base：folder
  # fileData:Data location
  #
  # Return：
  # "Data_prepro.csv"
  
  readSheet <- function(file){
    df <- read.csv(file,stringsAsFactors = F)
    rownames(df) <- df[,1]
    newdf <- df[,-1]
    return(newdf)
  }
  if (!all(file.exists(base))) 
    stop("'base' does not exists")
  info <- file.info(base)
  if (!all(info$isdir) && !all(!info$isdir)) 
    stop("'base needs to be either directories or files")
  if (all(info$isdir)) {
    csvfiles <- list.files(base, recursive = recursive, pattern = pattern, 
                           ignore.case = ignore.case, full.names = TRUE)
    if (verbose) {
      cat("[read] Found the following CSV files:\n")
      print(csvfiles)
    }
  }
  else csvfiles <- list.files(base, full.names = TRUE)
  dfs <- lapply(csvfiles, readSheet)
  namesUnion <- Reduce(union, lapply(dfs, names))
  df <- do.call(cbind, dfs)
  
  namedf <- names(df)
  dropname <- lapply(namedf, function(x){if(nchar(x) > 8) return(x)  
    else return(NULL)})
  dropname <- unlist(dropname)
  myname <- namedf %in% dropname
  df <- df[!myname]
  
  write.csv(df, file=paste(fileData, "Data_prepro.csv", sep = "/"))
  
  return(paste(fileData, "Data_prepro.csv", sep = "/"))
}

#base <- "E:/Rproject/ChAMP_data_DMRs/data_KIRC/afterComBatData"
FindExistCSV <- function(base, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = FALSE){
  # 函数说明：寻找在一个文件夹下存在的CSV文件，及批次名，afterComBatData_6042308144.csv
  #
  # Args：
  # base：位于那个文件夹下
  #
  # Return：
  # 6042308144, 6042308145
  
  GetBatchName <- function(filen){
    numstr <- str_extract_all(filen,"[0-9]+[0-9]")
    return(numstr[[1]])
  }
  if (!all(file.exists(base))) 
    stop("'base' does not exists")
  info <- file.info(base)
  if (!all(info$isdir) && !all(!info$isdir)) 
    stop("'base needs to be either directories or files")
  if (all(info$isdir)) {
    csvfiles <- list.files(base, recursive = recursive, pattern = pattern, 
                           ignore.case = ignore.case, full.names = TRUE)
    if (verbose) {
      cat("[read] Found the following CSV files:\n")
      print(csvfiles)
    }
  }
  else csvfiles <- list.files(base, full.names = TRUE)
  dfs <- lapply(csvfiles, GetBatchName)
  dfs <- unlist(dfs)
  return(as.numeric(dfs))
  
}



SureEachExamTarget <- function(eachExam){
  # 函数说明：Determine whether each sample is normal or cancer, cancer 1, normal -1
  #
  # Args：
  # eachExam：String type，eg. "T24"
  #
  # Returns：
  # target value：1 or -1
  
  if (substr(eachExam, 1, 1) == 'T') target <- 1
  else  target <- -1
  return(target)
}

ConvertPreproToModelData <- function(filePreproData,fileData){
  # Transform preprocessed data into data that can be recognized by the model
  #
  # Args：
  # filePreproData <- "resultsChamp/afterComBatData_bacth4.csv"
  # fileData: Data location
  #
  # Returns：
  # 'trainData.csv'，
  
  preData <- read.csv(filePreproData, stringsAsFactors =F)  
  rownames(preData) <- preData[,1]
  newData <- preData[,-1]
  rm(preData)
  newData.t <- t(newData)
  rm(newData)
  #newData <- as.data.frame(newData)
  exampleTN <- rownames(newData.t)
  target <- lapply(exampleTN,SureEachExamTarget)
  target <- unlist(target)
  newData.t <- data.frame(newData.t,target)
  newData <- t(newData.t)
  rm(newData.t)
  
  write.csv(newData,file = paste(fileData, 'trainData.csv', sep = "/"),row.names = T)
  fileTrain <- paste(fileData, 'trainData.csv', sep = "/")
  return(fileTrain)
}


UnionAndConvertDataToModel <- function(base, fileData, pattern = "csv$", ignore.case = TRUE, recursive = TRUE, verbose = TRUE){
  # Merge the data and transform it into model-recognizable data
  #
  # Args：
  # base：folder
  # fileData: Data location
  #
  # Return：
  # "Data_prepro.csv"
  
  readSheet <- function(file){
    df <- read.csv(file,stringsAsFactors = F)
    rownames(df) <- df[,1]
    newdf <- df[,-1]
    return(newdf)
  }
  if (!all(file.exists(base))) 
    stop("'base' does not exists")
  info <- file.info(base)
  if (!all(info$isdir) && !all(!info$isdir)) 
    stop("'base needs to be either directories or files")
  if (all(info$isdir)) {
    csvfiles <- list.files(base, recursive = recursive, pattern = pattern, 
                           ignore.case = ignore.case, full.names = TRUE)
    if (verbose) {
      cat("[read] Found the following CSV files:\n")
      print(csvfiles)
    }
  }
  else csvfiles <- list.files(base, full.names = TRUE)
  dfs <- lapply(csvfiles, readSheet)
  namesUnion <- Reduce(union, lapply(dfs, names))
  df <- do.call(cbind, dfs)
  
  namedf <- names(df)
  dropname <- lapply(namedf, function(x){if(nchar(x) > 8) return(x) 
    else return(NULL)})
  dropname <- unlist(dropname)
  myname <- namedf %in% dropname
  df <- df[!myname]
  
  write.csv(df, file=paste(fileData, "Data_prepro.csv", sep = "/"))
  
  newData.t <- t(df)
  rm(df)
  #newData <- as.data.frame(newData)
  exampleTN <- rownames(newData.t)
  target <- lapply(exampleTN,SureEachExamTarget)
  target <- unlist(target)
  newData.t <- data.frame(newData.t,target)
  trainData <- t(newData.t)
  rm(newData.t)
  
  save(trainData, file = paste(fileData, "trainData.Rdata", sep="/"))
  write.csv(trainData,file = paste(fileData, 'trainData.csv', sep = "/"),row.names = T)
  fileTrain <- paste(fileData, 'trainData.csv', sep = "/")
  return(fileTrain)
}
