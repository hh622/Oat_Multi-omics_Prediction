#===========================================================
#script02_DivPanel_multiomics_prediction
#HH 2021-07-31
#==========================================================

#------------------------------
# 2021-07-31
# models include:
# GBLUP, T, G+T, G+M, G+M+T

#clear workspace
rm(list=ls()) 

# set trait ID
i=1

#loading libraries
library(BGLR)
library(foreach)
library(doParallel)
source('../../BSFG/Estimate_gcor_prediction.R')
# install.packages("sommer")
# packageVersion("sommer")#4.1.0 
# note 4.1.1 does not work for my multi-kernel models

#set dirs
DIR_input <- "../01_Omicsdat/"
DIR_output <- "../02_OmisPred/"
DIR_output_SingleRun <- "../02_OmisPred/SingleRun/"

#if not exist, create dirs
if(!dir.exists(DIR_output)) dir.create(DIR_output)
if(!dir.exists(DIR_output_SingleRun)) dir.create(DIR_output_SingleRun)

#-------------------------------------------------------
#SECTION 1: load input data 
#-------------------------------------------------------
list.files(DIR_input)

# load pheno data
ls1=load(file = paste0(DIR_input,"DivPanel_pheno.RData"))
ls1 #"phedat"        "hsq_Div_pheno"

#load relationship matrix data
ls2=load(file = paste0(DIR_input,"DivPanel_Omicsdata_GRM_TRM_MRM_333x333.RData"))
ls2 #"Amat" "Dmat" "Emat" "TRM"  "MRM" 

#-------------------------------------------------
#SECTION 2: multiple omics prediction 
#-------------------------------------------------

#set parameters for loop
cycles=50 #nubmer of re-sampling
results <- NULL #collect results of cycles=100 runs of resampling
DT <- phedat #use which phenotypic dataset
kfold <- 2 #define K-fold cross-validation; 50:50 sampling

#set up parallel computation
Ncpu=40
cl <- makeCluster(Ncpu)
registerDoParallel(cl)

ptm <- proc.time() # Start the clock!

# r <- 1; nIter=200; burnIn=50 #testing
# for (i in 1: 3) {
for (i in 1: ncol(DT)){

  print(paste("###########runnning trait: ",colnames(DT)[i],"; colID=",i,"##########"))
  
  #get hsq for each trait
  hsq <- hsq_Div_pheno$hsq[hsq_Div_pheno$Traitname==colnames(DT)[i]]
  
  #prediction with BGLR for each trait
  results <-  foreach(r = 1:cycles, .combine=rbind, .inorder = TRUE) %dopar% {
    library(BGLR)
    nIter=20000
    burnIn=5000

    print(paste("###########runnning trait: ",colnames(DT)[i],"; colID=",i,"; runID=",r,"##########"))
    
    #-----------------------------------
    #step 0: set up traing and test sets
    #partition
    testSet <- sample(rownames(DT),round(nrow(DT)/kfold))
    trainSet <- setdiff(rownames(DT), testSet)
    
    #mask phenotypes for traning and test sets
    DT.tmp <- data.frame(trn=DT[,i], test=DT[,i])
    rownames(DT.tmp) <- rownames(DT)
    DT.tmp[testSet,"trn"] <- NA
    DT.tmp[!is.na(DT.tmp$trn),"test"] <- NA
    head(DT.tmp)
    # get observed phedat in test set excluding NAs
    y_test=DT.tmp[testSet,"test"]
    
    #--------------------------------------------------
    #step 1-1: single Kernel - G model(GBLUP)
    ETA=NULL; ETA=list(list(K=Amat,model='RKHS'))
    fmG=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
    
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmG$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpG = data.frame(TraitName=colnames(DT)[i],
                              runID=r,
                              model="G",Method = 'BSFG',
                              pearson = cor(fmG$yHat,DT.tmp$test, use = "complete")/sqrt(hsq),
                              g_cor  = cor(fmG$yHat,DT.tmp$test, use = "complete")/sqrt(hsq),
                              R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                              R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    #notes: for calculating R2a and R2b, make sure lines in y_test and yHat_test in the same order
    #i.e. identical(rownames(DT.tmp[testSet,]), rownames(yHat[testSet,,drop=F])) must be TRUE
    # calcualte R2 in testset 
    # R2a=cor(y,yhat,method='pearson')^2
    # R2b=1-sum((y-yHat)^2)/sum((y-mean(y))^2)
    # R2b actually compare y=yHat+e vs y=yBar
    # to see how much variance explained, should use R2a
    
    #--------------------------------------------------
    #step 1-2: single Kernel - T model
    ETA=NULL; ETA=list(list(K=TRM,model='RKHS'))
    fmT=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
     
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmT$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpT = data.frame(TraitName=colnames(DT)[i],
                              runID=r,
                              model="T",Method = 'BSFG',
                              pearson = cor(DT.tmp[testSet,"test"],yHat[testSet,], use ="complete")/sqrt(hsq),
                              g_cor = estimate_gcor(data.frame(ID=testSet, obs = DT.tmp[testSet,"test"],pred = yHat[testSet,]),
                                                    Knn=Amat[testSet,testSet],sKnn=NULL,method = 'MCMCglmm',normalize = T)[['g_cor']],
                              R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                              R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    plot(DT.tmp[testSet,"test"],yHat[testSet,])
    (sse_simpmod=sum((y_test-mean(y_test))^2))
    (sse_newmod=sum((fitted(lm(y_test~yHat_test))-y_test)^2))
    
    #--------------------------------------------------
    #step 1-3: single Kernel - M model
    ETA=NULL; ETA=list(list(K=MRM,model='RKHS'))
    fmM=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
     
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmM$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpM = data.frame(TraitName=colnames(DT)[i],
                              runID=r,
                              model="M",Method = 'BSFG',
                              pearson = cor(DT.tmp[testSet,"test"],yHat[testSet,], use ="complete")/sqrt(hsq),
                              g_cor = estimate_gcor(data.frame(ID=testSet,obs = DT.tmp[testSet,"test"],pred = yHat[testSet,]),
                                                    Knn=Amat[testSet,testSet],sKnn=NULL,method = 'MCMCglmm',normalize = T)[['g_cor']],
                              R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                              R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    
    #--------------------------------------------------
    #step 2-1: Two Kernel - A+D model
    ETA=NULL; ETA=list(list(K=Amat,model='RKHS'))
    ETA[[2]] <- list(K = Dmat, model = "RKHS")
    fmAD=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
    
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmAD$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpAD = data.frame(TraitName=colnames(DT)[i],
                               runID=r,
                               model="A+D",Method = 'BSFG',
                               pearson = cor(fmAD$yHat,DT.tmp$test, use = "complete")/sqrt(hsq),
                               g_cor  = cor(fmAD$yHat,DT.tmp$test, use = "complete")/sqrt(hsq),
                               R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                               R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    
    #--------------------------------------------------
    #step 2-2: Two Kernel - G+T model
    ETA=NULL; ETA=list(list(K=Amat,model='RKHS'))
    ETA[[2]] <- list(K = TRM, model = "RKHS")
    fmGT=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
    
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmGT$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpGT = data.frame(TraitName=colnames(DT)[i],
                               runID=r,
                               model="G+T",Method = 'BSFG',
                               pearson = cor(yHat[testSet,],DT.tmp[testSet,"test"], use="complete")/sqrt(hsq),
                               g_cor = estimate_gcor(data.frame(ID=testSet,obs = DT.tmp[testSet,"test"],pred = yHat[testSet,]),
                                                     Knn=Amat[testSet,testSet],sKnn=NULL,method = 'MCMCglmm',normalize = T)[['g_cor']],
                               R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                               R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    
    #--------------------------------------------------
    #step 2-3: Two Kernel - G+M model
    ETA=NULL; ETA=list(list(K=Amat,model='RKHS'))
    ETA[[2]] <- list(K = MRM, model = "RKHS")
    fmGM=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)

    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmGM$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpGM = data.frame(TraitName=colnames(DT)[i],
                               runID=r,
                               model="G+M",Method = 'BSFG',
                               pearson = cor(yHat[testSet,],DT.tmp[testSet,"test"], use="complete")/sqrt(hsq),
                               g_cor = estimate_gcor(data.frame(ID=testSet,obs = DT.tmp[testSet,"test"],pred = yHat[testSet,]),
                                                     Knn=Amat[testSet,testSet],sKnn=NULL,method = 'MCMCglmm',normalize = T)[['g_cor']],
                               R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                               R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    
    #--------------------------------------------------
    #step 3-1: Three Kernel - A+D+E model
    ETA=NULL; ETA=list(list(K=Amat,model='RKHS'))
    ETA[[2]] <- list(K = Dmat, model = "RKHS")
    ETA[[3]] <- list(K = Emat, model = "RKHS")
    fmADE=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
    
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmADE$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    results.tmpADE = data.frame(TraitName=colnames(DT)[i],
                                runID=r,
                                model="A+D+E",Method = 'BSFG',
                                pearson = cor(fmADE$yHat,DT.tmp$test, use = "complete")/sqrt(hsq),
                                g_cor  = cor(fmADE$yHat,DT.tmp$test, use = "complete")/sqrt(hsq),
                                R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                                R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    
    
    #--------------------------------------------------
    #step 3-2: Three Kernel - G+T+M model
    ETA=list(list(K=Amat,model='RKHS'))
    ETA[[2]] <- list(K = TRM, model = "RKHS")
    ETA[[3]] <- list(K = MRM, model = "RKHS")
    fmGTM=BGLR(y=DT.tmp$trn, ETA=ETA, nIter=nIter, burnIn=burnIn)
    
    #calculate prediction accuracy
    yHat=NULL; yHat <- data.frame(yHat=fmGTM$yHat,stringsAsFactors = F); rownames(yHat) <- rownames(DT.tmp)
    yHat_test=NULL; yHat_test=yHat[testSet,"yHat"]
    head(yHat)
    results.tmpGTM = data.frame(TraitName=colnames(DT)[i],
                                runID=r,
                                model="G+T+M",Method = 'BSFG',
                                pearson = cor(yHat[testSet,],DT.tmp[testSet,"test"], use="complete")/sqrt(hsq),
                                g_cor = estimate_gcor(data.frame(ID=testSet,obs = DT.tmp[testSet,"test"],pred = yHat[testSet,]),
                                                      Knn=Amat[testSet,testSet],sKnn=NULL,method = 'MCMCglmm',normalize = T)[['g_cor']],
                                R2a=summary(lm(y_test~yHat_test))$adj.r.squared,
                                R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2))
    
    #-------------------------------------------------------------
    #step 4: accumulate results from different runs of re-sampling
    results.tmp <- rbind(results.tmpG, results.tmpT, results.tmpM,
                         results.tmpAD, results.tmpGT, results.tmpGM,
                         results.tmpADE, results.tmpGTM)
    
    results.tmpG  <- NULL; results.tmpT<- NULL; results.tmpM<- NULL;
    results.tmpAD <- NULL; results.tmpGT <- NULL; results.tmpGM <- NULL;
    results.tmpADE <- NULL; results.tmpGTM <- NULL
    
    write.csv(results.tmp,
              paste0(DIR_output_SingleRun,
                     "DivPanel_multiomics_prediction_runID-",r,"_",colnames(DT)[i],".csv"))
    return(results.tmp)
    
    results.tmp <- NULL #empty results.tmp
  } #end of the foreach loop
  
  #output results
  write.csv(results,
            paste0(DIR_output,"DivPanel_multiomics_prediction_runs_",cycles,"_",colnames(DT)[i],".csv"))
}

stopCluster(cl)
proc.time() - ptm # Stop the clock
