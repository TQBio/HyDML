####################################################
# Programming instructions：Hybrid ensemble feature selection algorithm 
#                           for identifying differential methylation sites
# Final Edit Time ：2019.3.7
###################################################

# Data preprocessing removes sites with small variance
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


merge_my <- function(x,y){
  com <- merge(x,y,by="attribute")
  return(com)
}

#######################################################
OneTimeModel_FeatureSelectBySample_L2 <- function(x, preData, whichFS, underTimes.tum, underTimes.nor, 
                                                  overTimes,cv.fit){
  # From the data diversity training feature selection algorithm, find the differential methylation site
  #
  #
  #
  # return：
  # list(siteName=siteName, SummaryRSVM=fSummaryRSVM)
  time3 <- Sys.time()
  
  normalData <- subset(preData, target==-1)
  tumorData <- subset(preData, target==1)
  traintest <- UndersampledTumor(tumorData, normalData, underTimes.tum, underTimes.nor, overTimes)  
  testData_xy <- ConvertGlmnetData(traintest$testData)
  trainData_xy <- ConvertGlmnetData(traintest$trainData)
  feature_list <- list()
  acc_list <- list()
  ranking_list <- list()
  
  
  
  if("mcfs" == whichFS){
    #projections = "auto"
    model_mcfs <- mcfs(target~., traintest$trainData, featureFreq =100, projections = "auto",balance = 0,
                       cutoffPermutations = 0,threadsNumber=2)
    feature_mcfs <- subset(model_mcfs$RI,RI_norm>0,select=c(attribute,RI_norm))
    feature_mcfs <- feature_mcfs$attribute
    
    ranking_mcfs <- subset(model_mcfs$RI,select=c(attribute,RI_norm))
    ranking_mcfs$rank_mcfs <- ranking_mcfs$RI_norm
    ranking_mcfs <- ranking_mcfs[,-2]
    
    acc_mcfs <- ModelPerformance(feature_mcfs,testData_xy,trainData_xy)
    
    #     cat("mcfs:",'num_feature ：', length(feature_mcfs),'\n')
    #     cat("mcfs:","test_acc",acc_mcfs,"\n")
    
    feature_list <- feature_mcfs
    acc_list <- acc_mcfs
    ranking_list <- ranking_mcfs
  }
  
  if("glm" == whichFS){
    x <- trainData_xy$x
    y <- trainData_xy$y
    # cv.fit=cv.glmnet(x, y,family="binomial",alpha = 0.005)
    fit <- glmnet(x, y, family ="binomial",alpha = 0.005)
    coef_  <- coef(fit, s = cv.fit$lambda.min)
    
    #cat("--L2-lambda.min:",cv.fit$lambda.min,"\n")
    #lambda.min = cv.fit$lambda.min
    #cvm.min = min(cv.fit$cvm)
    #coef_ <- coef(cv.fit,s="lambda.min")
    rank_glm <- rep(0, dim(x)[2])
    index <- coef_@i[-1]
    rank_glm[index] <- abs(coef_@x[-1])
    attribute <- coef_@Dimnames[[1]][-1]
    
    feature_glm <- attribute[index] 
    acc_glm <- ModelPerformance(feature_glm,testData_xy,trainData_xy)
    rank_glm <- rank_glm*acc_glm[1]
    
    ranking_glm <- data.frame(attribute,rank_glm,stringsAsFactors = F)
    ranking_glm$attribute <- as.character(ranking_glm$attribute)
  
    feature_list <- feature_glm
    acc_list <- acc_glm
    ranking_list <- ranking_glm
  }
  

  if("MDFS" == whichFS){
    
    for (j in timesCross_L1) {
      
      feature_list[[j]] <- list()
      traintest <- Gettraintest(preData, underTimes.tum, underTimes.nor, overTimes)
      
      for (i in timesCross_L2) {
        
        
        cat("--L2--Time of find DML", i+j*10-10,"\n")
        
        trainData <- rbind(traintest[[1]],traintest[[2]],traintest[[3]],traintest[[4]],traintest[[5]])
        trainData <- trainData[-which(row.names(trainData)%in%row.names(traintest[[i]])),]
        testData <- traintest[[i]]
        
        traindata_xy <- ConvertGlmnetData(trainData)
        testdata_xy <- ConvertGlmnetData(testData)
        rm(trainData, testData)
        
        model_MDFS <- MDFS(traindata_xy$x,traindata_xy$y,dimensions = 2, divisions = 1)
        MDFS<-model_MDFS[[MDFS]]
        attribute<-row.names(preData) 
        rank<-MDFS[[IG]]
        ranking_MDFS<-data.frame(attribute,rank,stringsAsFactors = F)
        
        ranking_list[[(i+j*10-10)]] <- ranking_MDFS
        acc_MDFS <- ModelPerformance(ranking_MDFS,testdata_xy$x,traindata_xy$y)
        rank_MDFS <- rank_MDFS*acc_MDFS[1]
		feature_list <- attribute[index]
      }
      
    }
    feature <- feature_list
    ranking <- Merge_ranking(ranking_list)
    acc <- acc_list
  }
##############Function DIVERSITY
  
  times_results <- lapply(timesCross_L2, OneTimeModel_FeatureSelectBySample_L2, traintest$trainData,
                          whichFS, underTimes.tum,underTimes.nor,overTimes,cv.fit)
  
  ranking_list_L2 <- lapply(times_results, function(x){x$ranking})
  ranking <- Reduce( merge_my,ranking_list_L2)
  myvars <- colnames(ranking) %in% c("attribute")
  rank_ensemble <- apply(ranking[!myvars],1,sum) 
  #rank_ensemble <- apply(ranking[!myvars],1,function(x){length(x[which(x>0)])})  
  attribute <- ranking$attribute
  ranking_ensemble <- data.frame(attribute,rank_ensemble,stringsAsFactors = F)
  #feature_ensemble <- subset(ranking_ensemble, rank_ensemble>=mean(ranking_ensemble$rank_ensemble))$attribute
  feature_ensemble <- subset(ranking_ensemble, rank_ensemble>0)$attribute
  
  #cat("glm--ensemble:",'num_feature ：', length(feature_ensemble),'\n')
  acc_ensemble <- ModelPerformance(feature_ensemble,testData_xy,trainData_xy)
  #cat("glm--ensemble:","test_acc",acc_ensemble,"\n")
  
  
  feature_list$feature_ensemble <- feature_ensemble
  acc_list$acc_ensemble <- acc_ensemble
  ranking_list$ranking_ensemble <- ranking_ensemble
  
  ranking <- Reduce( merge_my,ranking_list)
  
  acc <- acc_list
  feature <- feature_list
  ############### aggregation
  
  time4 <- Sys.time()
  #cat("Primary modeling time：",time4-time3, "\n")
  
  return(list(ranking=ranking, acc=acc, feature=feature))
   
  

 
  return(list(ranking=ranking, acc=acc, feature=feature))
  
}

  if(flag_EFS){
    feature_EFS <- Reduce(intersect, feature_list)
    
    times_EFS <- table(unlist(feature_list))
    attribute <- colnames(trainData_xy$x)
    rank_EFS <- rep(0, dim(trainData_xy$x)[2])
    names(rank_EFS) <- attribute
    rank_EFS[names(times_EFS)] <- times_EFS
    
    ranking_EFS <- data.frame(attribute, rank_EFS,stringsAsFactors = F)
    acc_EFS <- ModelPerformance(feature_EFS,testData_xy,trainData_xy)
    cat("EFS:",'num_feature ：', length(feature_EFS),"\n")
    cat("EFS:","test_acc",acc_EFS,"\n")
    
    feature_list$feature_EFS <- feature_EFS
    acc_list$acc_EFS <- acc_EFS
    ranking_list$ranking_EFS <- ranking_EFS
  }
  

  
  ranking <- Reduce( merge_my,ranking_list)

  acc <- acc_list
  feature <- feature_list
  
  time4 <- Sys.time()
  cat("Time consuming：",time4-time3, "\n")
  
  return(list(ranking=ranking, acc=acc, feature=feature))
  
}

ModelPerformance <- function(feature,testData_xy,trainData_xy){
  # The features selected by the feature selection algorithm are evaluated by the test set
  # 
  # Args：
  # feature：
  # testData_xy：
  # trainData_xy ：
  
  svmres <- svm(trainData_xy$x[, feature], trainData_xy$y, type="C-classification",
                kernel="linear" )
  svmpred <- predict(svmres, testData_xy$x[, feature] )
  
  pred.y.frbs <- as.numeric(as.vector(svmpred))
  test.y.frbs <- as.numeric(as.vector(testData_xy$y))
  acc <- sum(pred.y.frbs == test.y.frbs )/dim(testData_xy$x)[1]  
  
#   evaluation index
#   tableOfFrbs <- table(as.matrix(pred.y.frbs), as.matrix(test.y.frbs))
#   #str(as.matrix(pred.y.frbs))
#   #str(as.matrix(test.y.frbs))
#   cat("tableOfFrbs:",tableOfFrbs,"\n")
#   TP <- tableOfFrbs[1, 1]
#   FP <- tableOfFrbs[1, 2]
#   FN <- tableOfFrbs[2, 1]
#   TN <- tableOfFrbs[2, 2]
#   ACCOfFrbs <- (TP + TN) / length(pred.y.frbs)
#   FPROfFrbs <-  FP / (TN + FP) 
#   TPROfFrbs <-  TP/(TP + FN)  
#   PrecisionFrbs <- TP/(TP + FP)
#   RecallFrbs <-  TP/(TP + FN)
#   F1Frbs <- 2*PrecisionFrbs*RecallFrbs/(PrecisionFrbs+RecallFrbs)
#   aucOfFrbs <- auc(as.vector(as.matrix(test.y.frbs)), as.vector(as.matrix(pred.y.frbs)))
#   
#   
#   performance_model <- c(ACCOfFrbs,FPROfFrbs,TPROfFrbs,PrecisionFrbs,RecallFrbs,F1Frbs,aucOfFrbs)
#   names(performance_model) <- c("ACC","FPR","TPR","Precision","Recall","F1","AUC")
  
  performance_model <- acc
  return(performance_model)
  
}
#########Sampling for DATA DIVERSITY
UndersampledTumor <- function(tumorData, normalData, underTimes.tum=0.6, underTimes.nor=0.6, overTimes=0){
  # Sampling of cancer samples and normal samples were combined into training samples, and the rest were test samples
  #
  #
  
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

ConvertGlmnetData <- function(dd ){
  # Convert the data to glmnet identifiable data
  # 
  # Args：
  # dd：ConvertPreproToModelData
  #
  # Returns：
  # list，x：feature ，y：label
  
  namesOfdd <- colnames(dd) %in% c('target')
  x <- as.matrix( dd[!namesOfdd] )
  target.dd <- dd$target
  #target.dd[which(target.dd==0)] = -1
  y <- factor( target.dd )
  return(list(x=x, y=y))
}


############# Feature selection robustness analysis #######################
GetRobustMetric <- function(oneFS, times_results_EFS){

  
  rank_oneFS <- paste0("rank_",oneFS)
  ranking_list <- lapply(times_results_EFS , function(x){x$ranking[,c("attribute",rank_oneFS)]})
  ranking <- Reduce( merge_my,ranking_list)
  ranking$attribute <- NULL
  
  jaccard_sub <- StabilityMetricOfFS(ranking, 'jaccard')
  robustMetric <- data.frame(jaccard_sub)
  rownames(robustMetric) <- oneFS
  
  cat(oneFS,"jaccard_sub:",jaccard_sub,"\n")
  
  return(robustMetric)
}

GetRobustMetric_meiguo <- function(oneFS, times_results_EFS){
  # The robustness of the feature selection algorithm is calculated


  
  rank_oneFS <- paste0("rank_",oneFS)
  ranking_list <- lapply(times_results_EFS , function(x){x$ranking[,c("attribute",rank_oneFS)]})
  ranking <- Reduce( merge_my,ranking_list)
  ranking$attribute <- NULL
  
  jaccard_sub <- StabilityMetricOfFS(ranking, 'jaccard')
  robustMetric <- data.frame(spearman_rank)
  rownames(robustMetric) <- oneFS
  
  cat(oneFS,"spearman_rank:",spearman_rank,"\n")
  
  return(robustMetric)
}

StabilityMetricOfFS <- function(rankingFeature, whichStabMetric){
  # Evaluation index of robustness of feature selection algorithm
  # 
  # Args：
  # rankingFeature: Feature importance ordering matrix, row represents the feature, and col represents the result of feature selection
  # whichStabMetric: jaccard,we use the Jaccard index
 
  
  change01 <- function(x){
    x[which(x>0)] <- 1
    return(x)
  }
  
  k <- dim(rankingFeature)[2]  
  comb <- combn(seq(1,k), 2) 
  num_comb <- k*(k-1)/2
  if(whichStabMetric == "jaccard"){
    subset <- apply(rankingFeature,2,change01)
    all_jaccard <- apply(comb, 2, StabilityMetric_jaccard, subset )
    stabilityMetric <- sum(all_jaccard)/num_comb
  }
  
  return(stabilityMetric)
  
}


StabilityMetric_jaccard <- function(oneComb, subset){
  # Calculate the similarity of two feature sets

  
  oneSub <- subset[,oneComb]
  oneSub <- t(oneSub)
  jaccard = 1 - distance(oneSub, method = "jaccard")
  
  return(jaccard)
  
}


###########Computational model performance #####################
GetModelMetric <- function(oneFS, times_results_EFS){
  # Analyze the performance of each feature selection algorithm model (cross-validation results)
  
  acc_oneFS <- paste0("acc_",oneFS)
  acc_list <- lapply(times_results_EFS , function(x){x$acc[[acc_oneFS]]})
  acc <- Reduce(rbind,acc_list)
  acc_mean <- apply(acc, 2, mean)
  modelMetric <- data.frame(acc_mean)
  colnames(modelMetric) <- oneFS
  
  cat(oneFS,acc_mean,"\n")
  
  return(modelMetric)
}



GetFeatureOfEFS <- function(oneFS="EFS", times_results_EFS, cutoff){
# Get the features obtained by the hybrid ensemble feature selection algorithm 
  
  feature_oneFS <- paste0("feature_",oneFS)
  feature_list <- lapply(times_results_EFS , function(x){x$feature[[feature_oneFS]]})
  
  feature <- table(unlist(feature_list))
  feature_f <- names(feature[which(feature>=cutoff)])
  
  return(feature_f)
}



######Unified naming & sorting######
rank_by_ranking<-function(ranking){
  ranking <- ranking[order(ranking[,2],decreasing=T),]
  ranking <- ranking[ranking_glm[,2]>0.001,]  
  return(ranking)
}

#######Normalization#############
sort_gyh<-function(x,ranking){
  a<-ranking$rank
  ranking[x,2]=(ranking[x,2]-ranking[length(ranking$attribute),2])/(ranking$rank[1]-ranking$rank[length(ranking$rank)])
  return(ranking)
}

########Take the union and rank it##
ranking_union<-function(ranking_MDFS,ranking_glm,ranking_caretRF){
  test1 <- intersect(ranking_MDFS$attribute,ranking_caretRF$attribute)
  test2 <- intersect(ranking_MDFS$attribute,ranking_glm$attribute)
  test3 <- intersect(ranking_glm$attribute,ranking_caretRF$attribute)
  a<-union(test1,test2)
  feature_name<-union(a,test3)
  
  data_feature<-data.frame(attribute=c(1:length(feature_name)),rank=c(1:length(feature_name)))
  data_feature$attribute<- feature_name
  # Rank merge
  for(i in 1:length(feature_name)){
    
    if(data_feature$attribute[i] %in% ranking_MDFS$attribute){
      t1=as.numeric(ranking_MDFS$rank[which(ranking_MDFS$attribute==data_feature$attribute[i])])
    } else {t1 =0}
    
    if(data_feature$attribute[i] %in% ranking_caretRF$attribute){
      t2=as.numeric(ranking_caretRF$rank[which(ranking_caretRF$attribute==data_feature$attribute[i])])
    }else {t2=0}
    
    if(data_feature$attribute[i] %in% ranking_glm$attribute){
      t3=as.numeric(ranking_glm$rank[which(ranking_glm$attribute==data_feature$attribute[i])])
    }else {t3=0}
    
    t=t1+t2+t3
    data_feature$rank[i]=as.numeric(t)
  }
  return(data_feature)
}

process<-function(ranking_MDFS,ranking_glm,ranking_caretRF){
  ranking_MDFS<-rank_by_ranking(ranking_MDFS)
  ranking_glm<-rank_by_ranking(ranking_glm)
  ranking_caretRF<-rank_by_ranking(ranking_caretRF)
  
  ranking_MDFS$rank<-lapply(1:length(ranking_MDFS$rank),sort_gyh,ranking_MDFS)
  ranking_caretRF$rank<-lapply(1:length(ranking_caretRF$rank),sort_gyh,ranking_caretRF)
  ranking_glm$rank<-lapply(1:length(ranking_glm$rank),sort_gyh,ranking_glm)
  
  data_fin<-ranking_union(ranking_MDFS,ranking_glm,ranking_caretRF)
  return(data_fin)
}


get_similarity<-function(BLCA_FIN,BRCA_FIN,COAD_FIN,
                         ESCA_FIN,HNSC_FIN,KIRC_FIN,
                         KIRP_FIN,LIHC_FIN,LUAD_FIN,
                         LUSC_FIN,PRAD_FIN,THCA_FIN,
                         UCEC_FIN){
  feature_list<-list(BLCA_FIN$attribute,BRCA_FIN$attribute,
                     COAD_FIN$attribute,ESCA_FIN$attribute,
                     HNSC_FIN$attribute,KIRC_FIN$attribute,
                     KIRP_FIN$attribute,LIHC_FIN$attribute,
                     LUAD_FIN$attribute,LUSC_FIN$attribute,
                     PRAD_FIN$attribute,THCA_FIN$attribute,
                     UCEC_FIN$attribute)
  
  feature<-table(unlist(feature_list))
  #fea_table <- sort(table(feature),decreasing=T)
  select_fea <- feature[which(feature>=10)]
  feature_similarity<-names(select_fea)
  write.table(feature,file ="feature.txt",sep=",")
  write.csv(feature_similarity, file = "similarity.csv")
  return(feature_similarity)
}




