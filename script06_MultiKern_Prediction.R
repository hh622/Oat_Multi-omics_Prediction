#==============================================
#script03_MultiKern_Prediction_in_Elite_Panel
#HH 2021-01-21
#==============================================

#-----------------
#2020-01-01
#predict each ENV seperately, 
#for each run, mask the same set of lines from each ENV

#-----------------
#2020-01-30
#based on network module eigengene and FA significant markers and markers in LD
#models: GBLUP, BayesB, MEs (three levels: Rsq>=0.1, 0.15, 0.20), FA models

library(BGLR)
library(foreach)
library(doSNOW)
library(reshape2)

DIR_input <- "../01_MultiKernPredInput/"
DIR_output <- "../02_output/ElitePanel/"
DIR_output_SingleRun <- "../02_output/ElitePanel/SingleRun/"
dir.create("../02_output/")
dir.create(DIR_output); dir.create(DIR_output_SingleRun)
dir.exists(DIR_output); dir.exists(DIR_output_SingleRun)

#-------------------------------------------------------------------------
# grab the array id value from the environment variable passed from sbatch
#-------------------------------------------------------------------------
# slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
# 
# # coerce the value to an integer
# i <- as.numeric(slurm_arrayid)
# # i=10

#if run each trait separately on gore02
args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1]) #for testing script
# i = 3
#-------------------------------
# Load phedat and GRM

lsdata=load(paste0(DIR_input,"MultiKernPred_GRMs_ElitePanel.RData"))
lsdata
# [1] "phedat2"             "genodat"             "trait.H2.acrosssite" "trait.H2.singlesite"
# [5] "G_all"               "G01_MEs_Rsq0.10"     "G02_MEs_Rsq0.10"     "G01_MEs_Rsq0.15"    
# [9] "G02_MEs_Rsq0.15"     "G01_MEs_Rsq0.20"     "G02_MEs_Rsq0.20"     "G01_FA"             
# [13] "G02_FA" 

dim(G_all)#211 211
dim(genodat)#211 71682

#extract only FA traits
phedat=phedat2 #rename phedat2 as phedat
trait.H2 <- trait.H2.acrosssite
head(phedat)
identical(trait.H2$Traitname, 
          colnames(phedat)[-which(colnames(phedat2) %in%c("LOC","LINE"))]) #TRUE
phedat <- phedat[,c(1:2,which(grepl("FAME.",colnames(phedat))))]
trait.H2 <- trait.H2[grepl("FAME.",trait.H2$Traitname),]
identical(colnames(phedat)[-1*1:2],trait.H2$Traitname) #TRUE
dim(phedat) #633  12
head(phedat)
which(grepl("FAME.",colnames(phedat)))
# [1]  3  4  5  6  7  8  9 10 11 12

#--------------------------------------------------------------------------------
#prepare phenotypes - dcast one vector to three columns (each column is one ENV)
Y0 <- phedat[,c(1,2,i)]
Y <- dcast(Y0, LINE ~ LOC, value.var = colnames(phedat)[i])
rownames(Y) <- Y$LINE
Y <- Y[,-which(colnames(Y)=="LINE")]
Y <- as.matrix(Y)
head(Y)
str(Y)
# 
#-------------------------------
#set up running pars
rundate <- "20210130"
kfold <- 5
ncycles=50 #nubmer of re-sampling
nIter=20000; burnIn=5000
DT0 <- Y #phedat without masking used for validation later on
DT = Y  #phedat will be masked for model fitting 
results <- NULL
results.tmp <- NULL
head(DT)

#create computational cluster manually
cl <- makeCluster(ncycles, type = "SOCK")
registerDoSNOW(cl)
# cl <- makeCluster(ncycles)
# registerDoParallel(cl)

#get hsq used for calculating prediction accuracy
hsq <- trait.H2$hsq[trait.H2$Traitname==colnames(phedat)[i]]

#get n and nEnv
n <- nrow(Y);  nEnv <- ncol(Y)

head(DT)
# r <- 1
ptm <- proc.time() # Start the clock!
results =  foreach(r = 1:ncycles, .combine=rbind, .inorder = TRUE, .errorhandling="pass") %dopar% {
  library(BGLR)
  library(reshape2)
  # nIter=20000; burnIn=5000
  nIter=200; burnIn=50 #for testing
  
  #set up Env names
  env1 <- "MN"; env2 <- "SD"; env3 <- "WI"
  
  print(paste("###########runnning trait: ",colnames(DT)[i],"; colID=",i,"; runID=",r,"##########"))
  
  #----------------------------------
  #step 0: set up trn and test sets
  #partition
  testSet <- sample(rownames(DT),round(nrow(DT)/kfold))
  trainSet <- setdiff(rownames(DT), testSet)
  length(testSet); length(trainSet) #42; 169
  DT[testSet,] <- NA
  colSums(is.na(DT))
  # MN SD WI 
  # 42 45 47 
  
  #-------------------------------------
  #model 1: GBLUP
  eigen_G_all <- eigen(G_all)
  ETA <- NULL; ETA <- list(G=list(V=eigen_G_all$vectors,d=eigen_G_all$values,model='RKHS'))
  fm_GBLUPEnv1 <-BGLR(y=DT[,env1],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_GBLUPEnv2 <-BGLR(y=DT[,env2],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_GBLUPEnv3 <-BGLR(y=DT[,env3],ETA=ETA,nIter=nIter,burnIn=burnIn)
  Predability_GBLUPEnv1 = cor(DT0[testSet,env1],fm_GBLUPEnv1$yHat[testSet],use = "na.or.complete")
  Predability_GBLUPEnv2 = cor(DT0[testSet,env2],fm_GBLUPEnv2$yHat[testSet],use = "na.or.complete")
  Predability_GBLUPEnv3 = cor(DT0[testSet,env3],fm_GBLUPEnv3$yHat[testSet],use = "na.or.complete")
  
  #-------------------------------------
  #model 2: BayesB
  ETA <- NULL; ETA <- list(G=list(X=genodat, model='BayesB'))
  fm_BayesBEnv1 <-BGLR(y=DT[,env1],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_BayesBEnv2 <-BGLR(y=DT[,env2],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_BayesBEnv3 <-BGLR(y=DT[,env3],ETA=ETA,nIter=nIter,burnIn=burnIn)
  Predability_BayesBEnv1 = cor(DT0[testSet,env1],fm_BayesBEnv1$yHat[testSet],use = "na.or.complete")
  Predability_BayesBEnv2 = cor(DT0[testSet,env2],fm_BayesBEnv2$yHat[testSet],use = "na.or.complete")
  Predability_BayesBEnv3 = cor(DT0[testSet,env3],fm_BayesBEnv3$yHat[testSet],use = "na.or.complete")
  
  #------------------------------------------------------------
  #model 3: MEs_Rsq0.10 ([1] "G01_MEs_Rsq0.10" "G02_MEs_Rsq0.10" )
  eigenG01_MEs_Rsq0.10 <- eigen(G01_MEs_Rsq0.10)
  eigenG02_MEs_Rsq0.10 <- eigen(G02_MEs_Rsq0.10)
  
  ETA <- NULL
  ETA <- list(G=list(V=eigenG01_MEs_Rsq0.10$vectors,d=eigenG01_MEs_Rsq0.10$values,model='RKHS'))
  ETA[[2]] <- list(V=eigenG02_MEs_Rsq0.10$vectors,d=eigenG02_MEs_Rsq0.10$values,model='RKHS')
  fm_MEs_Rsq0.10Env1 <-BGLR(y=DT[,env1],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MEs_Rsq0.10Env2 <-BGLR(y=DT[,env2],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MEs_Rsq0.10Env3 <-BGLR(y=DT[,env3],ETA=ETA,nIter=nIter,burnIn=burnIn)
  Predability_MEs_Rsq0.10Env1 = cor(DT0[testSet,env1],fm_MEs_Rsq0.10Env1$yHat[testSet],use = "na.or.complete")
  Predability_MEs_Rsq0.10Env2 = cor(DT0[testSet,env2],fm_MEs_Rsq0.10Env2$yHat[testSet],use = "na.or.complete")
  Predability_MEs_Rsq0.10Env3 = cor(DT0[testSet,env3],fm_MEs_Rsq0.10Env3$yHat[testSet],use = "na.or.complete")
  
  #-----------------------------------------------------------
  #model 4: MEs_Rsq0.15 ([1] "G01_MEs_Rsq0.15" "G02_MEs_Rsq0.15")
  eigenG01_MEs_Rsq0.15 <- eigen(G01_MEs_Rsq0.15)
  eigenG02_MEs_Rsq0.15 <- eigen(G02_MEs_Rsq0.15)
  ETA <- NULL
  ETA <- list(G=list(V=eigenG01_MEs_Rsq0.15$vectors,d=eigenG01_MEs_Rsq0.15$values,model='RKHS'))
  ETA[[2]] <- list(V=eigenG02_MEs_Rsq0.15$vectors,d=eigenG02_MEs_Rsq0.15$values,model='RKHS')
  fm_MEs_Rsq0.15Env1 <-BGLR(y=DT[,env1],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MEs_Rsq0.15Env2 <-BGLR(y=DT[,env2],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MEs_Rsq0.15Env3 <-BGLR(y=DT[,env3],ETA=ETA,nIter=nIter,burnIn=burnIn)
  Predability_MEs_Rsq0.15Env1 = cor(DT0[testSet,env1],fm_MEs_Rsq0.15Env1$yHat[testSet],use = "na.or.complete")
  Predability_MEs_Rsq0.15Env2 = cor(DT0[testSet,env2],fm_MEs_Rsq0.15Env2$yHat[testSet],use = "na.or.complete")
  Predability_MEs_Rsq0.15Env3 = cor(DT0[testSet,env3],fm_MEs_Rsq0.15Env3$yHat[testSet],use = "na.or.complete")
  
  #------------------------------------------------------------
  #model 5: MEs_Rsq0.20 ([1] "G01_MEs_Rsq0.20" "G02_MEs_Rsq0.20")
  eigenG01_MEs_Rsq0.20 <- eigen(G01_MEs_Rsq0.20)
  eigenG02_MEs_Rsq0.20 <- eigen(G02_MEs_Rsq0.20)
  
  ETA <- NULL
  ETA <- list(G=list(V=eigenG01_MEs_Rsq0.20$vectors,d=eigenG01_MEs_Rsq0.20$values,model='RKHS'))
  ETA[[2]] <- list(V=eigenG02_MEs_Rsq0.20$vectors,d=eigenG02_MEs_Rsq0.20$values,model='RKHS')
  fm_MEs_Rsq0.20Env1 <-BGLR(y=DT[,env1],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MEs_Rsq0.20Env2 <-BGLR(y=DT[,env2],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MEs_Rsq0.20Env3 <-BGLR(y=DT[,env3],ETA=ETA,nIter=nIter,burnIn=burnIn)
  Predability_MEs_Rsq0.20Env1 = cor(DT0[testSet,env1],fm_MEs_Rsq0.20Env1$yHat[testSet],use = "na.or.complete")
  Predability_MEs_Rsq0.20Env2 = cor(DT0[testSet,env2],fm_MEs_Rsq0.20Env2$yHat[testSet],use = "na.or.complete")
  Predability_MEs_Rsq0.20Env3 = cor(DT0[testSet,env3],fm_MEs_Rsq0.20Env3$yHat[testSet],use = "na.or.complete")
  
  #-------------------------------------
  #model 6: FA facotor informative model
  eigenG01_MKlip <- eigen(G01_FA)
  eigenG02_MKlip <- eigen(G02_FA)
  
  ETA <- NULL
  ETA <- list(G=list(V=eigenG01_MKlip$vectors,d=eigenG01_MKlip$values,model='RKHS'))
  ETA[[2]] <- list(V=eigenG02_MKlip$vectors,d=eigenG02_MKlip$values,model='RKHS')
  fm_MKlipEnv1 <-BGLR(y=DT[,env1],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MKlipEnv2 <-BGLR(y=DT[,env2],ETA=ETA,nIter=nIter,burnIn=burnIn)
  fm_MKlipEnv3 <-BGLR(y=DT[,env3],ETA=ETA,nIter=nIter,burnIn=burnIn)
  Predability_MKlipEnv1 = cor(DT0[testSet,env1],fm_MKlipEnv1$yHat[testSet],use = "na.or.complete")
  Predability_MKlipEnv2 = cor(DT0[testSet,env2],fm_MKlipEnv2$yHat[testSet],use = "na.or.complete")
  Predability_MKlipEnv3 = cor(DT0[testSet,env3],fm_MKlipEnv3$yHat[testSet],use = "na.or.complete")
  
  
  #---------------------------------------
  # rbind prediction data
  results.tmp <- data.frame(Kernel="G",
                            Fold=kfold,
                            TraitName=colnames(phedat)[i],
                            Model = c(rep("GBLUP",3),
                                      rep("BayesB",3),
                                      rep("MEs_Rsq0.10",3),
                                      rep("MEs_Rsq0.15",3),
                                      rep("MEs_Rsq0.20",3), 
                                      rep("MKlip",3)),
                            TestENV=rep(c(env1,env2,env3),6),
                            RunID=r,
                            PredAbility=c(Predability_GBLUPEnv1,
                                          Predability_GBLUPEnv2,
                                          Predability_GBLUPEnv3,
                                          Predability_BayesBEnv1,
                                          Predability_BayesBEnv2,
                                          Predability_BayesBEnv3,
                                          Predability_MEs_Rsq0.10Env1,
                                          Predability_MEs_Rsq0.10Env2,
                                          Predability_MEs_Rsq0.10Env3,
                                          Predability_MEs_Rsq0.15Env1,
                                          Predability_MEs_Rsq0.15Env2,
                                          Predability_MEs_Rsq0.15Env3,
                                          Predability_MEs_Rsq0.20Env1,
                                          Predability_MEs_Rsq0.20Env2,
                                          Predability_MEs_Rsq0.20Env3,
                                          Predability_MKlipEnv1,
                                          Predability_MKlipEnv2,
                                          Predability_MKlipEnv3),
                            PredAccuracy=c(Predability_GBLUPEnv1,
                                           Predability_GBLUPEnv2,
                                           Predability_GBLUPEnv3,
                                           Predability_BayesBEnv1,
                                           Predability_BayesBEnv2,
                                           Predability_BayesBEnv3,
                                           Predability_MEs_Rsq0.10Env1,
                                           Predability_MEs_Rsq0.10Env2,
                                           Predability_MEs_Rsq0.10Env3,
                                           Predability_MEs_Rsq0.15Env1,
                                           Predability_MEs_Rsq0.15Env2,
                                           Predability_MEs_Rsq0.15Env3,
                                           Predability_MEs_Rsq0.20Env1,
                                           Predability_MEs_Rsq0.20Env2,
                                           Predability_MEs_Rsq0.20Env3,
                                           Predability_MKlipEnv1,
                                           Predability_MKlipEnv2,
                                           Predability_MKlipEnv3)/sqrt(hsq),
                            stringsAsFactors = F)
  #output results of a single run
  write.csv(results.tmp,
            paste0(DIR_output_SingleRun,"MultiKern_Prediction_ElitePanel_TaskID-",
                   i,"_",colnames(phedat)[i],"_",kfold,"fold_","rundate",rundate,"_run-",r,".csv"),
            row.names = F)
  
  return(results.tmp)
  results.tmp <- NULL #empty results.tmp
  
}#end of foreach loop parallel run

#output results
write.csv(results,
          paste0(DIR_output,"MultiKern_Prediction_ElitePanel_TaskID-", i,"_runs_",ncycles,"_",
                 colnames(phedat)[i],".csv"))

stopCluster(cl)
proc.time() - ptm # Stop the clock


