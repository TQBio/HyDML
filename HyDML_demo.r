library(ecodist)  
library(e1071)
library(glmnet)
library(rmcfs)
library(MDFS)
library(sampling)
library(pROC)
library(randomForest)

####parameter####
fileData <- "C:/Users/asus/Desktop/R/data"  # data file location and you can download training data from 
                                            # https://github.com/TQBio/HyDML/tree/master/Data
filew <- fileData
underTimes.tum <- 0.9
underTimes.nor <- 0.9 
overTimes <- 0  
timesCross_L1 <- 1:5 
timesCross_L2 <- 1:5

####subfunction####
RemoveSmallVar <- function(preData, pct=1.0/5){
  # The sites with small variance after preprocessing were filtered out
  # 
  # Args：
  # preData：data
  # pct：Percentage of filtered data defult: 0.2
  #
  # Return:
  # varPreData：The data with smaller variance points were filtered out
  
  var.data <- apply(preData,1,var)
  #var.data <- unlist(var.data)
  
  newvar <- var.data[!(names(var.data) %in% c("target"))]
  newvar <- sort(newvar)
  
  dropSite <-  round(length(newvar) * pct)
  nameNeedSite <- names(newvar[dropSite:length(newvar)])
  nameNeedSite <- c(nameNeedSite, "target")
  varPreData <- preData[nameNeedSite,]
  return(varPreData)
  
}

UndersampledTumor <- function(tumorData, normalData, underTimes.tum=0.9, underTimes.nor=0.9,
                              overTimes=0){
  # After undersampling, the cancer samples were combined with normal samples to form training samples, and the rest were test samples
  #
  # Args：
  # tumorData：cancer samples
  # normalData：normal sample
  # underTimes.tum: the proportion of undersampled cancer samples in the training samples to the total cancer samples
  # underTimes.nor: the proportion of the undersampled normal samples in the training samples to the total normal samples
  # overTimes: oversampling normal sample : training cancer samples = 1 : overtimes. if overtimes == 0, no oversampling
  #
  # return：
  # list(testData=testData, trainData=trainData)
  
  len.tumor <- dim(tumorData)[1]
  len.norma <- dim(normalData)[1]
  
  SampInd <- seq(1, len.tumor)
  TrainInd <- sample(SampInd, floor(underTimes.tum*len.tumor), replace=F )
  TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
  
  
  SampInd.norL1 <- seq(1, len.norma)
  TrainInd.norL1 <- sample(SampInd.norL1, floor(underTimes.nor*len.norma), replace=F )
  TestInd.norL1 <- SampInd.norL1[ which(!(SampInd.norL1 %in% TrainInd.norL1 ))]
  
  normalData.train <- normalData[TrainInd.norL1, ]
  
  lentrainTum <- length(TrainInd)
  nOver <- floor(lentrainTum*overTimes) - length(TrainInd.norL1)
  if(nOver > 0){
    
    TrainInd.nor <- sample(TrainInd.norL1, nOver, replace=T )  
    #TestInd <- SampInd[ which(!(SampInd %in% TrainInd ))]
    
    #testData <- normalData[TestInd,]
    normalData.new <- rbind(normalData.train, normalData[TrainInd.nor,])
    normalData.train <- normalData.new
    rm(normalData.new)
    
  }
  
  
  testData <- rbind(tumorData[TestInd,], normalData[TestInd.norL1, ])
  trainData <- rbind(tumorData[TrainInd,], normalData.train)
  
  #   cat("train_tumor:", length(TrainInd), "train_normal:",dim(normalData.train)[1], "\n")
  #   cat("test_tumor:", length(TestInd), "test_normal:",length(TestInd.norL1),"\n")
  #   
  
  return(list(testData=testData, trainData=trainData))
}

Gettraintest <- function(preData, underTimes.tum, underTimes.nor, overTimes){
  # preData was randomly divided into training samples and test samples
  #
  # Args：
  # preData：data
  # underTimes.tum: the proportion of undersampled cancer samples in the training samples to the total cancer samples
  # underTimes.nor: the proportion of the undersampled normal samples in the training samples to the total normal samples
  # overTimes: oversampling normal sample : training cancer samples = 1 : overtimes. if overtimes == 0, no oversampling
  #
  # return：
  # traintest：a set of training and test samples for each group
  
  normalData <- subset(preData, target==-1)
  a <- floor(length(normalData[,1])/5)
  normaldata <- list()
  normaldata <- split.data.frame(normalData, sample(rep(1:5, c(a,a,a,a,a))))
  rm(normalData)
  
  tumorData <- subset(preData, target==1)
  b <- floor(length(tumorData[,1])/5)
  tumordata <- list()
  tumordata <- split.data.frame(tumorData, sample(rep(1:5, c(b,b,b,b,b))))
  rm(tumorData)
  
  traintest <- list()
  for(i in 1:5){
    
    traintest[[i]] <- UndersampledTumor(tumordata[[i]], normaldata[[i]], underTimes.tum, underTimes.nor, overTimes)
    traintest[[i]] <- rbind(traintest[[i]][["trainData"]], traintest[[i]][["testData"]])
    
  }
  
  return(traintest)
  
}


ConvertGlmnetData <- function(dd ){
  # convert the data to glmnet and MDFS identifiable data
  # 
  # Args：
  # dd：data
  #
  # Returns：
  # list(x=x, y=y)
  # x: feature set
  # y: label
  
  namesOfdd <- colnames(dd) %in% c('target')
  x <- as.matrix( dd[!namesOfdd] )
  target.dd <- dd$target
  #target.dd[which(target.dd==0)] = -1
  y <- factor( target.dd )
  return(list(x=x, y=y))
}

ModelPerformance <- function(traindata_xy, testdata_xy, classifier, feature){
  # The features selected by the feature selection algorithm are evaluated by the test set
  # 
  # Args：
  # testData_xy：test set
  # trainData_xy ：train set
  # classifier：select "svm" or "randomforest" to classify
  # feature：the differential methylation site
  #
  # return：
  # auc: area Under Curve
  # acc: accuracy
  
  if(classifier=="svm"){
    
    svmres <- svm(traindata_xy$x[, feature], traindata_xy$y, type="C-classification",
                  kernel="linear", cost=100, gamma=1)
    svmpred <- predict(svmres, testdata_xy$x[, feature] )
    
    pred.y.frbs <- as.numeric(as.vector(svmpred))
    test.y.frbs <- as.numeric(as.vector(testdata_xy$y))
    
    #      tableOfFrbs <- table(as.matrix(pred.y.frbs), as.matrix(test.y.frbs))
    #      TP <- tableOfFrbs[1, 1]
    #      FP <- tableOfFrbs[1, 2]
    #      FN <- tableOfFrbs[2, 1]
    #      TN <- tableOfFrbs[2, 2]
    #      SPOfFrbs <- (1 - FP / (TN + FP))  
    #      SEOfFrbs <-  TP/(TP + FN) 
    auc <- auc(as.vector(as.matrix(test.y.frbs)), as.vector(as.matrix(pred.y.frbs)))
    #      acc <- sum(pred.y.frbs == test.y.frbs )/dim(testdata_xy$x)[1]
    
  }
  
  if(classifier=="randomforest"){
    
    rf_model<-randomForest(traindata_xy$x[, feature], traindata_xy$y,
                           importance = TRUE, proximity = TRUE)
    rf_pre<-predict(rf_model,testdata_xy$x[, feature])
    
    pred.y.frbs <- as.numeric(as.vector(rf_pre))
    test.y.frbs <- as.numeric(as.vector(testdata_xy$y))
    
    #      tableOfFrbs <- table(as.matrix(pred.y.frbs), as.matrix(test.y.frbs))
    #      TP <- tableOfFrbs[1, 1]
    #      FP <- tableOfFrbs[1, 2]
    #      FN <- tableOfFrbs[2, 1]
    #      TN <- tableOfFrbs[2, 2]
    #      SPOfFrbs <- (1 - FP / (TN + FP))  
    #      SEOfFrbs <-  TP/(TP + FN)  
    auc <- auc(as.vector(as.matrix(test.y.frbs)), as.vector(as.matrix(pred.y.frbs)))
    #      acc  <- sum(pred.y.frbs == test.y.frbs )/dim(testdata_xy$x)[1] 
    
  }
  
  return(auc) # or "acc"
  
}

FeatureSelect_glm <- function(preData, timesCross_L1, timesCross_L2,
                              underTimes.tum, underTimes.nor, overTimes){
  # Using "glmnet" for feature selection algorithm to feature select
  #
  # Args:
  # preData：data
  # timesCross_L1, timesCross_L2: the times of cross-validation
  # underTimes.tum: the proportion of undersampled cancer samples in the training samples to the total cancer samples
  # underTimes.nor: the proportion of the undersampled normal samples in the training samples to the total normal samples
  # overTimes: oversampling normal sample : training cancer samples = 1 : overtimes. if overtimes == 0, no oversampling
  # 
  # Return:
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM, auc=auc)
  
  feature_list <- list()
  ranking_list <- list()
  auc <- list()
  for (j in timesCross_L1) {
    
    feature_list[[j]] <- list()
    traintest <- Gettraintest(preData, underTimes.tum, underTimes.nor, overTimes)

    for (i in timesCross_L2) {

      cat("times: ", i+j*5-5,"\n")
      
      trainData <- rbind(traintest[[1]],traintest[[2]],traintest[[3]],traintest[[4]],traintest[[5]])
      trainData <- trainData[-which(row.names(trainData)%in%row.names(traintest[[i]])),]
      testData <- traintest[[i]]
      traindata_xy <- ConvertGlmnetData(trainData)
      testdata_xy <- ConvertGlmnetData(testData)
      rm(trainData)
      
      cv.fit=cv.glmnet(traindata_xy$x, traindata_xy$y,family="binomial",alpha = 0.005)
      fit <- glmnet(traindata_xy$x, traindata_xy$y, family ="binomial",alpha = 0.005)
      coef_  <- coef(fit, s = cv.fit$lambda.min)
      rank_glm <- rep(0, dim(traindata_xy$x)[2])
      index <- coef_@i[-1]
      rank_glm[index] <- abs(coef_@x[-1])
      attribute <- coef_@Dimnames[[1]][-1]
      
      feature_glm <- attribute[index] 
      ranking_glm <- data.frame(attribute,rank,stringsAsFactors = F)
      ranking_glm$attribute <- as.character(ranking_glm$attribute)
      
      feature_list[[j]][[i]] <- feature_glm
      ranking_list[[(i+j*5-5)]] <- ranking_glm
      auc[[(i+j*5-5)]] <- ModelPerformance(traindata_xy, testdata_xy, "svm", feature_glm)
      
    }
    
  }
  return(list(feature_list, ranking_list, auc))
}

FeatureSelect_MDFS <- function(preData, timesCross_L1, timesCross_L2,
                              underTimes.tum, underTimes.nor, overTimes){
  # Using "MDFS" for feature selection algorithm to feature select
  #
  # Args:
  # preData：data
  # timesCross_L1, timesCross_L2: the times of cross-validation
  # underTimes.tum: the proportion of undersampled cancer samples in the training samples to the total cancer samples
  # underTimes.nor: the proportion of the undersampled normal samples in the training samples to the total normal samples
  # overTimes: oversampling normal sample : training cancer samples = 1 : overtimes. if overtimes == 0, no oversampling
  # 
  # Return:
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM)
  
  feature_list <- list()
  ranking_list <- list()
  auc <- list()
  for (j in timesCross_L1) {
    
    feature_list[[j]] <- list()
    traintest <- Gettraintest(preData, underTimes.tum, underTimes.nor, overTimes)
    
    for (i in timesCross_L2) {
      
      cat("times: ", i+j*5-5,"\n")
      
      trainData <- rbind(traintest[[1]],traintest[[2]],traintest[[3]],traintest[[4]],traintest[[5]])
      trainData <- trainData[-which(row.names(trainData)%in%row.names(traintest[[i]])),]
      testData <- traintest[[i]]
      traindata_xy <- ConvertGlmnetData(trainData)
      testdata_xy <- ConvertGlmnetData(testData)
      rm(trainData)
      
      model_MDFS <- MDFS(traindata_xy$x,traindata_xy$y,dimensions = 2, divisions = 1)
      MDFS<-model_MDFS[[MDFS]]
      attribute<-row.names(preData) 
      rank<-MDFS[[IG]]
      ranking_MDFS<-data.frame(attribute,rank,stringsAsFactors = F)
      feature_MDFS <- subset(ranking_MDFS,rank_MDFS>0)$attribute
      
      
      feature_list[[j]][[i]] <- feature_MDFS
      ranking_list[[(i+j*5-5)]] <- ranking_MDFS
      auc[[(i+j*5-5)]] <- ModelPerformance(traindata_xy, testdata_xy, "svm", feature_MDFS)
      
    }
    
  }
  return(list(feature_list, ranking_list, auc))
}

FeatureSelect_rmcfs <- function(preData, timesCross_L1, timesCross_L2,
                              underTimes.tum, underTimes.nor, overTimes){
  # Using "mcfs" for feature selection algorithm to feature select
  #
  # Args:
  # preData：data
  # timesCross_L1, timesCross_L2: the times of cross-validation
  # underTimes.tum: the proportion of undersampled cancer samples in the training samples to the total cancer samples
  # underTimes.nor: the proportion of the undersampled normal samples in the training samples to the total normal samples
  # overTimes: oversampling normal sample : training cancer samples = 1 : overtimes. if overtimes == 0, no oversampling
  # 
  # Return:
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM)
  
  feature_list <- list()
  ranking_list <- list()
  auc <- list()
  for (j in timesCross_L1) {
    
    feature_list[[j]] <- list()
    traintest <- Gettraintest(preData, underTimes.tum, underTimes.nor, overTimes)
    
    for (i in timesCross_L2) {
      
      cat("times: ", i+j*5-5,"\n")
      
      trainData <- rbind(traintest[[1]],traintest[[2]],traintest[[3]],traintest[[4]],traintest[[5]])
      trainData <- trainData[-which(row.names(trainData)%in%row.names(traintest[[i]])),]
      testData <- traintest[[i]]
      
      model_mcfs <- mcfs(target~., trainData, featureFreq =100, projections = "auto",balance = 0,
                         cutoffPermutations = 0,threadsNumber=2)
      feature_mcfs <- subset(model_mcfs$RI,RI_norm>0,select=c(attribute,RI_norm))
      feature_mcfs <- feature_mcfs$attribute
      
      ranking_mcfs <- subset(model_mcfs$RI,select=c(attribute,RI_norm))
      ranking_mcfs$rank_mcfs <- ranking_mcfs$RI_norm
      ranking_mcfs <- ranking_mcfs[,-2]
      
      feature_list[[j]][[i]] <- feature_mcfs
      ranking_list[[(i+j*5-5)]] <- ranking_mcfs
      auc[[(i+j*5-5)]] <- ModelPerformance(traindata_xy, testdata_xy, "svm", feature_mcfs)
      
    }
    
  }
  return(list(feature_list, ranking_list, auc))
}

Merge_ranking <- function(ranking_list){
  # get the union of all the rank in the ranking_list
  # 
  # Args：
  # ranking_list：the set of feature ranking
  #
  # return:
  # ranking: feature ranking
  
  ranking <- merge(ranking_list[[1]], ranking_list[[2]], by="attribute",all = TRUE)
  ranking <- data.frame(attribute=ranking$attribute, rank=ranking$rank.x + ranking$rank.y)
  
  for (k in 3:25) {
    
    ranking <- merge(ranking, ranking_list[[k]], by="attribute",all = TRUE)
    ranking <- data.frame(attribute=ranking$attribute, rank=ranking$rank + ranking$rank)
    
  }
  
  return(ranking)
  
}

rank_by_ranking<-function(ranking){
  # sort features by ranking
  # 
  # Args：
  # ranking: feature ranking
  #
  # Return:
  # ranking: feature ranking
  
  ranking <- ranking[order(ranking[,2],decreasing=T),]
  ranking <- ranking[ranking_glm[,2]>0.001,]  
  return(ranking)
}

sort_gyh<-function(x,ranking){
  # normalize the ranking values obtained by three methods
  # 
  # Args：
  # x: the number of ranking
  # ranking: feature ranking
  #
  # Return:
  # ranking: feature ranking
  
  a<-ranking$rank
  ranking[x,2]=(ranking[x,2]-ranking[length(ranking$attribute),2])/(ranking$rank[1]-ranking$rank[length(ranking$rank)])
  return(ranking)
  
}

ranking_union<-function(ranking_glm,ranking_mcfs,ranking_MDFS){
  # combine the ranking values obtained by three methods
  # 
  # Args：
  # ranking_glm: feature ranking by using "glm"
  # ranking_mcfs: feature ranking by using "MDFS"
  # ranking_MDFS: feature ranking by using "mcfs"
  #
  # Return:
  # data_feature: the differential methylation site
  
  test1 <- intersect(ranking_MDFS$attribute,ranking_mcfs$attribute)
  test2 <- intersect(ranking_MDFS$attribute,ranking_glm$attribute)
  test3 <- intersect(ranking_glm$attribute,ranking_mcfs$attribute)
  a<-union(test1,test2)
  feature_name<-union(a,test3)
  
  data_feature<-data.frame(attribute=c(1:length(feature_name)),rank=c(1:length(feature_name)))
  data_feature$attribute<- feature_name

    for(i in 1:length(feature_name)){
    
    if(data_feature$attribute[i] %in% ranking_MDFS$attribute){
      t1=as.numeric(ranking_MDFS$rank[which(ranking_MDFS$attribute==data_feature$attribute[i])])
    } else {t1 =0}
    
    if(data_feature$attribute[i] %in% ranking_mcfs$attribute){
      t2=as.numeric(ranking_mcfs$rank[which(ranking_mcfs$attribute==data_feature$attribute[i])])
    }else {t2=0}
    
    if(data_feature$attribute[i] %in% ranking_glm$attribute){
      t3=as.numeric(ranking_glm$rank[which(ranking_glm$attribute==data_feature$attribute[i])])
    }else {t3=0}
    
    t=t1+t2+t3
    data_feature$rank[i]=as.numeric(t)
  }
  return(data_feature)
}

process<-function(ranking_glm,ranking_mcfs,ranking_MDFSF){
  # the general process of findingthe differential methylation sites
  # 
  # Args：
  # ranking_glm: feature ranking by using "glm"
  # ranking_mcfs: feature ranking by using "MDFS"
  # ranking_MDFS: feature ranking by using "mcfs"
  #
  # Return:
  # data_fin: the differential methylation site
  
  ranking_glm<-rank_by_ranking(ranking_glm)
  ranking_mcfs<-rank_by_ranking(ranking_mcfs)
  ranking_MDFS<-rank_by_ranking(ranking_MDFS)
  
  ranking_glm$rank<-lapply(1:length(ranking_glm$rank),sort_gyh,ranking_glm)
  ranking_mcfs$rank<-lapply(1:length(ranking_mcfs$rank),sort_gyh,ranking_mcfs)
  ranking_MDFS$rank<-lapply(1:length(ranking_MDFS$rank),sort_gyh,ranking_MDFS)
  
  data_fin<-ranking_union(ranking_glm,ranking_mcfs,ranking_MDFSF)
  return(data_fin)
}

####main function####
start <- Sys.time()

fileTrainData <- paste(fileData,'trainData.Rdata', sep = "/")
load(fileTrainData)

cat("before var select data length: ", dim(trainData), '\n')
preData.dropSmallVar <- RemoveSmallVar(trainData, pct=1.0/2)  
rm(trainData)
preData <- t(preData.dropSmallVar)
preData <- as.data.frame(preData)
preData$target <- as.factor(preData$target)
rm(preData.dropSmallVar)
cat("after var select data length: ", dim(preData), '\n')

result_glm <- FeatureSelect_glm(preData, timesCross_L1, timesCross_L2,
                                underTimes.tum, underTimes.nor, overTimes)
result_mcfs <- FeatureSelect_mcfs(preData, timesCross_L1, timesCross_L2,
                                    underTimes.tum, underTimes.nor, overTimes)
result_MDFS <- FeatureSelect_MDFS(preData, timesCross_L1, timesCross_L2,
                                  underTimes.tum, underTimes.nor, overTimes)

ranking_glm <- Merge_ranking(result_glm[[2]])
ranking_mcfs <- Merge_ranking(result_mcfs[[2]])
ranking_MDFS <- Merge_ranking(result_MDFS[[2]])

data_fin <- process(ranking_glm,ranking_mcfs,ranking_MDFSF)

save(data_fin,file = file.path(filew,"data_fin.Rdata"))

end <- Sys.time()
cat("run_time:",end-start,'\n')
