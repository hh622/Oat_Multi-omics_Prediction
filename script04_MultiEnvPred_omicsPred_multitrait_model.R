#====================================================
#script02_ElitePanel_omicsPred_multitrait_model
#HH 2021-08-01
#====================================================

#HH 2021-08-01
# ElitePanel_omicsPred_multitrait_model
# models: 
# 1) SNP-based GBLUP
# 2) metabolites-based MBLUP
# 3) SNPs + Metabolites 2-kernel model
# reference: https://github.com/gdlc/BGLR-R (PART 2. Multi-trait models)

library(doParallel)
library(foreach)
library(BGLR)
library(dplyr)

#set dirs
DIR_input <- "../01_Omicsdat/"
DIR_output <- "../02_OmisPred/"
DIR_output_SingleRun <- "../02_OmisPred/SingleRun/"

#if not exist, create dirs
if(!dir.exists(DIR_output)) dir.create(DIR_output)
if(!dir.exists(DIR_output_SingleRun)) dir.create(DIR_output_SingleRun)

#-------------------------
#SECTION 01: read in data
#-------------------------
ls=load(file = paste0(DIR_input,"ElitePanel_OmicsPred_Input_GC02LC02.RData"))
ls
# [1] "phedat"              "trait.H2.acrosssite" "trait.H2.singlesite" "Ksnp"               
# [5] "KM_MN_SD"            "KM_MN_WI"            "KM_SD_WI"  
sum(is.na(Ksnp)) #0
sum(is.na(KM_MN_SD)) #0
sum(is.na(KM_MN_WI)) #0
sum(is.na(KM_MN_WI)) #0

#-------------------------------------
#SECTION 02: MET-GBLUP prediction
#-------------------------------------

#set up parallel computation
Ncpu=40
cl <- makeCluster(Ncpu)
registerDoParallel(cl)

#set up running pars
rundate <- "20210801"
kfold <- 5
cycles=50 #nubmer of re-sampling
pheno.startID <- grep("Plantheight",colnames(phedat)) 
pheno.endID <- length(colnames(phedat))
PhenoID_lst <- pheno.startID:pheno.endID
runstart <- 1; runend <- cycles
nIter=20000; burnIn=5000
results.tmp <- NULL
results <- NULL
results.all <- NULL

ptm <- proc.time() # Start the clock!

#run prediction with 50 re-samplings
test.script = FALSE
if(test.script){
  i<- pheno.startID; r <- 1; runstart =1; runend=2; nIter=100; burnIn=20 #testing
}
# for (i in 1: 3) {
for(i in pheno.startID:pheno.endID){
  print(paste("PhenoID =",i,"is started!")) #computation started HERE!!!
  
  #-------------------------------
  #prepare phenotypes
  head(phedat)
  DT.MN <- phedat[phedat$LOC=="MN",c("LINE",colnames(phedat)[i])]
  colnames(DT.MN) <- c("LINE","MN")
  head(DT.MN); dim(DT.MN)
  
  DT.SD <- phedat[phedat$LOC=="SD",c("LINE",colnames(phedat)[i])]
  colnames(DT.SD) <- c("LINE","SD")
  head(DT.SD); dim(DT.SD)
  
  DT.WI <- phedat[phedat$LOC=="WI",c("LINE",colnames(phedat)[i])]
  colnames(DT.WI) <- c("LINE","WI")
  head(DT.WI); dim(DT.WI)
  
  DT.oat <- left_join(left_join(DT.MN, DT.SD),DT.WI)
  head(DT.oat)
  
  #get hsq for trait i
  hsq <- trait.H2.acrosssite$hsq[trait.H2.acrosssite$Traitname==colnames(phedat)[i]]
  
  #------------------------------------
  #prediction with BGLR for each trait
  results=foreach(r = runstart:runend, .combine=rbind, .inorder = TRUE, .errorhandling="pass")  %dopar% {
  # for(r in runstart:runend) {
    print(paste("#######TaskID =",i,", foreach loop is started!######"))
    
    library(BGLR)
    library(reshape2)
    nIter=20000; burnIn=5000
    
    #----------------------------------
    #prepare Y.trn and Y.test of CV2
    #I will mask round(nrow(Y)/kfold)*3 missing data points
    #with the restriction that one genotype is missing in only one environment
    
    Y <- DT.oat
    rownames(Y) <- Y$LINE
    Y <- Y[,c("MN","SD","WI")]
    Y <- Y[rownames(Ksnp),] #re-order rows of Y according to LINE in GRM
    dim(Y)
    head(Y)
    cor(Y,use = "na.or.complete" )
    
    #double check phedat, GRM and MRM all in the same order
    sapply(list(rownames(Ksnp),colnames(Ksnp),
                rownames(KM_MN_SD),colnames(KM_MN_SD),
                rownames(KM_MN_WI),colnames(KM_MN_WI),
                rownames(KM_SD_WI),colnames(KM_SD_WI)), 
           FUN = identical, 
           rownames(Y)) #all returned TRUE
    
    #mask 1/Kfold phenotypic values as missings for each ENV separately
    tst <- sample(rownames(Y),round(nrow(Y)/kfold))
    Y.trn.MN <- as.matrix(Y)
    Y.trn.SD <- as.matrix(Y)
    Y.trn.WI <- as.matrix(Y)
    Y.trn.MN[rownames(Y.trn.MN) %in% tst, "MN"] <- NA
    Y.trn.SD[rownames(Y.trn.SD) %in% tst, "SD"] <- NA
    Y.trn.WI[rownames(Y.trn.WI) %in% tst, "WI"] <- NA
    
    #-----------------------------------------
    # fitting G-by-E models with G Kernel
    
    #1) DIAG-DIAG model
    ETA <- NULL
    ETA <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")))
    Gkern_MN_D_D <- Multitrait(y=Y.trn.MN, ETA=ETA,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    Gkern_SD_D_D <- Multitrait(y=Y.trn.SD, ETA=ETA,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    Gkern_WI_D_D <- Multitrait(y=Y.trn.WI, ETA=ETA,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    
    #2) DIAG-UN model
    ETA <- NULL
    ETA <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")))
    Gkern_MN_D_UN <- Multitrait(y=Y.trn.MN, ETA=ETA, nIter=nIter,burnIn=burnIn)#model residuals are UNstructured by default
    Gkern_SD_D_UN <- Multitrait(y=Y.trn.SD, ETA=ETA, nIter=nIter,burnIn=burnIn)
    Gkern_WI_D_UN <- Multitrait(y=Y.trn.WI, ETA=ETA, nIter=nIter,burnIn=burnIn)
    
    #3) UN-DIAG model
    ETA <- NULL
    ETA <-list(list(K=Ksnp,model="RKHS"))#the random effect is UNstructured by default
    Gkern_MN_UN_D <-Multitrait(y=Y.trn.MN, ETA=ETA,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Gkern_SD_UN_D <-Multitrait(y=Y.trn.SD, ETA=ETA,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Gkern_WI_UN_D <-Multitrait(y=Y.trn.WI, ETA=ETA,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    
    #4) UN-UN model (# the default setting of BGLR Multitrait: UN-UN)
    # the covariance matrix of the random effect and that of model residuals are UNstructured by default.
    ETA <- NULL
    ETA <-list(list(K=Ksnp,model="RKHS"))
    Gkern_MN_UN_UN <-Multitrait(y=Y.trn.MN, ETA=ETA,nIter=nIter,burnIn=burnIn)
    Gkern_SD_UN_UN <-Multitrait(y=Y.trn.SD, ETA=ETA,nIter=nIter,burnIn=burnIn)
    Gkern_WI_UN_UN <-Multitrait(y=Y.trn.WI, ETA=ETA,nIter=nIter,burnIn=burnIn)
    
    #5) FA-DIAG model
    M <- matrix(nrow = 3, ncol = 1, TRUE)
    ETA <- NULL
    ETA <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)))
    Gkern_MN_FA_D <-Multitrait(y=Y.trn.MN, ETA=ETA,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Gkern_SD_FA_D <-Multitrait(y=Y.trn.SD, ETA=ETA,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Gkern_WI_FA_D <-Multitrait(y=Y.trn.WI, ETA=ETA,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    
    #6) FA-UN model
    M <- matrix(nrow = 3, ncol = 1, TRUE)
    ETA <- NULL
    ETA <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)))
    Gkern_MN_FA_UN <- Multitrait(y=Y.trn.MN, ETA=ETA, nIter=nIter,burnIn=burnIn)#the random effect is UNstructured by default
    Gkern_SD_FA_UN <- Multitrait(y=Y.trn.SD, ETA=ETA, nIter=nIter,burnIn=burnIn)
    Gkern_WI_FA_UN <- Multitrait(y=Y.trn.WI, ETA=ETA, nIter=nIter,burnIn=burnIn)
    
    #-----------------------------------------
    # fitting G-by-E models with M Kernel
    #1) DIAG-DIAG model
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=KM_SD_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA2 <-list(list(K=KM_MN_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA3 <-list(list(K=KM_MN_SD,model="RKHS",Cov=list(type="DIAG")))
    Mkern_MN_D_D <- Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    Mkern_SD_D_D <- Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    Mkern_WI_D_D <- Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    
    #2) DIAG-UN model
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=KM_SD_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA2 <-list(list(K=KM_MN_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA3 <-list(list(K=KM_MN_SD,model="RKHS",Cov=list(type="DIAG")))
    Mkern_MN_D_UN <- Multitrait(y=Y.trn.MN, ETA=ETA1, nIter=nIter,burnIn=burnIn)
    Mkern_SD_D_UN <- Multitrait(y=Y.trn.SD, ETA=ETA2, nIter=nIter,burnIn=burnIn)
    Mkern_WI_D_UN <- Multitrait(y=Y.trn.WI, ETA=ETA3, nIter=nIter,burnIn=burnIn)
    
    #3) UN-DIAG model
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=KM_SD_WI,model="RKHS"))
    ETA2 <-list(list(K=KM_MN_WI,model="RKHS"))
    ETA3 <-list(list(K=KM_MN_SD,model="RKHS"))
    Mkern_MN_UN_D <-Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Mkern_SD_UN_D <-Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Mkern_WI_UN_D <-Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    
    #4) UN-UN model (# the default setting of BGLR Multitrait: UN-UN)
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=KM_SD_WI,model="RKHS"))
    ETA2 <-list(list(K=KM_MN_WI,model="RKHS"))
    ETA3 <-list(list(K=KM_MN_SD,model="RKHS"))
    Mkern_MN_UN_UN <-Multitrait(y=Y.trn.MN, ETA=ETA1,nIter=nIter,burnIn=burnIn)
    Mkern_SD_UN_UN <-Multitrait(y=Y.trn.SD, ETA=ETA2,nIter=nIter,burnIn=burnIn)
    Mkern_WI_UN_UN <-Multitrait(y=Y.trn.WI, ETA=ETA3,nIter=nIter,burnIn=burnIn)
    
    #5) FA-DIAG model
    M <- matrix(nrow = 3, ncol = 1, TRUE)
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=KM_SD_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA2 <-list(list(K=KM_MN_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA3 <-list(list(K=KM_MN_SD, model="RKHS", Cov=list(type="FA",M=M)))
    Mkern_MN_FA_D <-Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Mkern_SD_FA_D <-Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    Mkern_WI_FA_D <-Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    
    #6) FA-UN model
    M <- matrix(nrow = 3, ncol = 1, TRUE)
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=KM_SD_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA2 <-list(list(K=KM_MN_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA3 <-list(list(K=KM_MN_SD, model="RKHS", Cov=list(type="FA",M=M)))
    Mkern_MN_FA_UN <- Multitrait(y=Y.trn.MN, ETA=ETA1, nIter=nIter,burnIn=burnIn)
    Mkern_SD_FA_UN <- Multitrait(y=Y.trn.SD, ETA=ETA2, nIter=nIter,burnIn=burnIn)
    Mkern_WI_FA_UN <- Multitrait(y=Y.trn.WI, ETA=ETA3, nIter=nIter,burnIn=burnIn)
    
    #-----------------------------------------
    # fitting G-by-E models with GM Kernel
    
    #1) DIAG-DIAG model
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")),list(K=KM_SD_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA2 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")),list(K=KM_MN_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA3 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")),list(K=KM_MN_SD,model="RKHS",Cov=list(type="DIAG")))
    GMkern_MN_D_D <- Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    GMkern_SD_D_D <- Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    GMkern_WI_D_D <- Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="DIAG"),nIter=nIter,burnIn=burnIn)
    
    #2) DIAG-UN model
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")),list(K=KM_SD_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA2 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")),list(K=KM_MN_WI,model="RKHS",Cov=list(type="DIAG")))
    ETA3 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="DIAG")),list(K=KM_MN_SD,model="RKHS",Cov=list(type="DIAG")))
    GMkern_MN_D_UN <- Multitrait(y=Y.trn.MN, ETA=ETA1, nIter=nIter,burnIn=burnIn)
    GMkern_SD_D_UN <- Multitrait(y=Y.trn.SD, ETA=ETA2, nIter=nIter,burnIn=burnIn)
    GMkern_WI_D_UN <- Multitrait(y=Y.trn.WI, ETA=ETA3, nIter=nIter,burnIn=burnIn)
    
    #3) UN-DIAG model
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=Ksnp,model="RKHS"),list(K=KM_SD_WI,model="RKHS"))
    ETA2 <-list(list(K=Ksnp,model="RKHS"),list(K=KM_MN_WI,model="RKHS"))
    ETA3 <-list(list(K=Ksnp,model="RKHS"),list(K=KM_MN_SD,model="RKHS"))
    GMkern_MN_UN_D <-Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    GMkern_SD_UN_D <-Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    GMkern_WI_UN_D <-Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    
    #4) UN-UN model (# the default setting of BGLR Multitrait: UN-UN)
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=Ksnp,model="RKHS"),list(K=KM_SD_WI,model="RKHS"))
    ETA2 <-list(list(K=Ksnp,model="RKHS"),list(K=KM_MN_WI,model="RKHS"))
    ETA3 <-list(list(K=Ksnp,model="RKHS"),list(K=KM_MN_SD,model="RKHS"))
    GMkern_MN_UN_UN <-Multitrait(y=Y.trn.MN, ETA=ETA1,nIter=nIter,burnIn=burnIn)
    GMkern_SD_UN_UN <-Multitrait(y=Y.trn.SD, ETA=ETA2,nIter=nIter,burnIn=burnIn)
    GMkern_WI_UN_UN <-Multitrait(y=Y.trn.WI, ETA=ETA3,nIter=nIter,burnIn=burnIn)
    
    #5) FA-DIAG model
    M <- matrix(nrow = 3, ncol = 1, TRUE)
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)),list(K=KM_SD_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA2 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)),list(K=KM_MN_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA3 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)),list(K=KM_MN_SD, model="RKHS", Cov=list(type="FA",M=M)))
    GMkern_MN_FA_D <-Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    GMkern_SD_FA_D <-Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    GMkern_WI_FA_D <-Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="DIAG"), nIter=nIter,burnIn=burnIn)
    
    #6) FA-UN model
    M <- matrix(nrow = 3, ncol = 1, TRUE)
    ETA1 <- NULL; ETA2 <- NULL; ETA3 <- NULL
    ETA1 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)),list(K=KM_SD_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA2 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)),list(K=KM_MN_WI, model="RKHS", Cov=list(type="FA",M=M)))
    ETA3 <-list(list(K=Ksnp,model="RKHS",Cov=list(type="FA",M=M)),list(K=KM_MN_SD, model="RKHS", Cov=list(type="FA",M=M)))
    GMkern_MN_FA_UN <- Multitrait(y=Y.trn.MN, ETA=ETA1,resCov=list(type="UN"), nIter=nIter,burnIn=burnIn)
    GMkern_SD_FA_UN <- Multitrait(y=Y.trn.SD, ETA=ETA2,resCov=list(type="UN"), nIter=nIter,burnIn=burnIn)
    GMkern_WI_FA_UN <- Multitrait(y=Y.trn.WI, ETA=ETA3,resCov=list(type="UN"), nIter=nIter,burnIn=burnIn)
    
    #-------------------------------------------------------
    #calculate prediction ability only for masked data points
    #G kernel
    test.env1 <- "MN"; test.env2 <- "SD"; test.env3 <- "WI"
    Predability_Gkern_MN_D_D <- cor(Gkern_MN_D_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Gkern_MN_D_UN <- cor(Gkern_MN_D_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Gkern_MN_UN_D <- cor(Gkern_MN_UN_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Gkern_MN_UN_UN <- cor(Gkern_MN_UN_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Gkern_MN_FA_D <- cor(Gkern_MN_FA_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Gkern_MN_FA_UN <- cor(Gkern_MN_FA_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    
    Predability_Gkern_SD_D_D <- cor(Gkern_SD_D_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Gkern_SD_D_UN <- cor(Gkern_SD_D_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Gkern_SD_UN_D <- cor(Gkern_SD_UN_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Gkern_SD_UN_UN <- cor(Gkern_SD_UN_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Gkern_SD_FA_D <- cor(Gkern_SD_FA_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Gkern_SD_FA_UN <- cor(Gkern_SD_FA_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    
    Predability_Gkern_WI_D_D <- cor(Gkern_WI_D_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Gkern_WI_D_UN <- cor(Gkern_WI_D_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Gkern_WI_UN_D <- cor(Gkern_WI_UN_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Gkern_WI_UN_UN <- cor(Gkern_WI_UN_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Gkern_WI_FA_D <- cor(Gkern_WI_FA_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Gkern_WI_FA_UN <- cor(Gkern_WI_FA_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    
    #M kernel
    Predability_Mkern_MN_D_D <- cor(Mkern_MN_D_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Mkern_MN_D_UN <- cor(Mkern_MN_D_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Mkern_MN_UN_D <- cor(Mkern_MN_UN_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Mkern_MN_UN_UN <- cor(Mkern_MN_UN_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Mkern_MN_FA_D <- cor(Mkern_MN_FA_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_Mkern_MN_FA_UN <- cor(Mkern_MN_FA_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    
    Predability_Mkern_SD_D_D <- cor(Mkern_SD_D_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Mkern_SD_D_UN <- cor(Mkern_SD_D_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Mkern_SD_UN_D <- cor(Mkern_SD_UN_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Mkern_SD_UN_UN <- cor(Mkern_SD_UN_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Mkern_SD_FA_D <- cor(Mkern_SD_FA_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_Mkern_SD_FA_UN <- cor(Mkern_SD_FA_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    
    Predability_Mkern_WI_D_D <- cor(Mkern_WI_D_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Mkern_WI_D_UN <- cor(Mkern_WI_D_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Mkern_WI_UN_D <- cor(Mkern_WI_UN_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Mkern_WI_UN_UN <- cor(Mkern_WI_UN_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Mkern_WI_FA_D <- cor(Mkern_WI_FA_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_Mkern_WI_FA_UN <- cor(Mkern_WI_FA_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    
    #GM kernel
    Predability_GMkern_MN_D_D <- cor(GMkern_MN_D_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_GMkern_MN_D_UN <- cor(GMkern_MN_D_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_GMkern_MN_UN_D <- cor(GMkern_MN_UN_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_GMkern_MN_UN_UN <- cor(GMkern_MN_UN_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_GMkern_MN_FA_D <- cor(GMkern_MN_FA_D$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    Predability_GMkern_MN_FA_UN <- cor(GMkern_MN_FA_UN$yHat[tst,test.env1], Y[tst,test.env1], use = "na.or.complete")
    
    Predability_GMkern_SD_D_D <- cor(GMkern_SD_D_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_GMkern_SD_D_UN <- cor(GMkern_SD_D_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_GMkern_SD_UN_D <- cor(GMkern_SD_UN_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_GMkern_SD_UN_UN <- cor(GMkern_SD_UN_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_GMkern_SD_FA_D <- cor(GMkern_SD_FA_D$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    Predability_GMkern_SD_FA_UN <- cor(GMkern_SD_FA_UN$yHat[tst,test.env2], Y[tst,test.env2], use = "na.or.complete")
    
    Predability_GMkern_WI_D_D <- cor(GMkern_WI_D_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_GMkern_WI_D_UN <- cor(GMkern_WI_D_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_GMkern_WI_UN_D <- cor(GMkern_WI_UN_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_GMkern_WI_UN_UN <- cor(GMkern_WI_UN_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_GMkern_WI_FA_D <- cor(GMkern_WI_FA_D$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    Predability_GMkern_WI_FA_UN <- cor(GMkern_WI_FA_UN$yHat[tst,test.env3], Y[tst,test.env3], use = "na.or.complete")
    
    #----------------------------------------------------------------
    #calculate R2a (R2a=summary(lm(y_test~yHat_test))$adj.r.squared))$adj.r.squared$adj.r.squared))$adj.r.squared
    
    #G kernel
    R2a_Gkern_MN_D_D <- summary(lm(Gkern_MN_D_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Gkern_MN_D_UN <- summary(lm(Gkern_MN_D_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Gkern_MN_UN_D <- summary(lm(Gkern_MN_UN_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Gkern_MN_UN_UN <- summary(lm(Gkern_MN_UN_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Gkern_MN_FA_D <- summary(lm(Gkern_MN_FA_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Gkern_MN_FA_UN <- summary(lm(Gkern_MN_FA_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    
    R2a_Gkern_SD_D_D <- summary(lm(Gkern_SD_D_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Gkern_SD_D_UN <- summary(lm(Gkern_SD_D_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Gkern_SD_UN_D <- summary(lm(Gkern_SD_UN_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Gkern_SD_UN_UN <- summary(lm(Gkern_SD_UN_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Gkern_SD_FA_D <- summary(lm(Gkern_SD_FA_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Gkern_SD_FA_UN <- summary(lm(Gkern_SD_FA_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    
    R2a_Gkern_WI_D_D <- summary(lm(Gkern_WI_D_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Gkern_WI_D_UN <- summary(lm(Gkern_WI_D_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Gkern_WI_UN_D <- summary(lm(Gkern_WI_UN_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Gkern_WI_UN_UN <- summary(lm(Gkern_WI_UN_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Gkern_WI_FA_D <- summary(lm(Gkern_WI_FA_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Gkern_WI_FA_UN <- summary(lm(Gkern_WI_FA_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    
    #M kernel
    R2a_Mkern_MN_D_D <- summary(lm(Mkern_MN_D_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Mkern_MN_D_UN <- summary(lm(Mkern_MN_D_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Mkern_MN_UN_D <- summary(lm(Mkern_MN_UN_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Mkern_MN_UN_UN <- summary(lm(Mkern_MN_UN_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Mkern_MN_FA_D <- summary(lm(Mkern_MN_FA_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_Mkern_MN_FA_UN <- summary(lm(Mkern_MN_FA_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    
    R2a_Mkern_SD_D_D <- summary(lm(Mkern_SD_D_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Mkern_SD_D_UN <- summary(lm(Mkern_SD_D_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Mkern_SD_UN_D <- summary(lm(Mkern_SD_UN_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Mkern_SD_UN_UN <- summary(lm(Mkern_SD_UN_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Mkern_SD_FA_D <- summary(lm(Mkern_SD_FA_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_Mkern_SD_FA_UN <- summary(lm(Mkern_SD_FA_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    
    R2a_Mkern_WI_D_D <- summary(lm(Mkern_WI_D_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Mkern_WI_D_UN <- summary(lm(Mkern_WI_D_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Mkern_WI_UN_D <- summary(lm(Mkern_WI_UN_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Mkern_WI_UN_UN <- summary(lm(Mkern_WI_UN_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Mkern_WI_FA_D <- summary(lm(Mkern_WI_FA_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_Mkern_WI_FA_UN <- summary(lm(Mkern_WI_FA_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    
    #GM kernel
    R2a_GMkern_MN_D_D <- summary(lm(GMkern_MN_D_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_GMkern_MN_D_UN <- summary(lm(GMkern_MN_D_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_GMkern_MN_UN_D <- summary(lm(GMkern_MN_UN_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_GMkern_MN_UN_UN <- summary(lm(GMkern_MN_UN_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_GMkern_MN_FA_D <- summary(lm(GMkern_MN_FA_D$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    R2a_GMkern_MN_FA_UN <- summary(lm(GMkern_MN_FA_UN$yHat[tst,test.env1]~ Y[tst,test.env1]))$adj.r.squared
    
    R2a_GMkern_SD_D_D <- summary(lm(GMkern_SD_D_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_GMkern_SD_D_UN <- summary(lm(GMkern_SD_D_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_GMkern_SD_UN_D <- summary(lm(GMkern_SD_UN_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_GMkern_SD_UN_UN <- summary(lm(GMkern_SD_UN_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_GMkern_SD_FA_D <- summary(lm(GMkern_SD_FA_D$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    R2a_GMkern_SD_FA_UN <- summary(lm(GMkern_SD_FA_UN$yHat[tst,test.env2]~ Y[tst,test.env2]))$adj.r.squared
    
    R2a_GMkern_WI_D_D <- summary(lm(GMkern_WI_D_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_GMkern_WI_D_UN <- summary(lm(GMkern_WI_D_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_GMkern_WI_UN_D <- summary(lm(GMkern_WI_UN_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_GMkern_WI_UN_UN <- summary(lm(GMkern_WI_UN_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_GMkern_WI_FA_D <- summary(lm(GMkern_WI_FA_D$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    R2a_GMkern_WI_FA_UN <- summary(lm(GMkern_WI_FA_UN$yHat[tst,test.env3]~ Y[tst,test.env3]))$adj.r.squared
    
    #---------------------------------------------------------------------------
    #calculate R2b; R2b=1-sum((y_test-yHat_test)^2)/sum((y_test-mean(y_test))^2)
    
    #G kernel
    R2b_Gkern_MN_D_D <- 1-sum((Gkern_MN_D_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Gkern_MN_D_UN <- 1-sum((Gkern_MN_D_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Gkern_MN_UN_D <- 1-sum((Gkern_MN_UN_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Gkern_MN_UN_UN <- 1-sum((Gkern_MN_UN_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Gkern_MN_FA_D <- 1-sum((Gkern_MN_FA_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Gkern_MN_FA_UN <- 1-sum((Gkern_MN_FA_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    
    R2b_Gkern_SD_D_D <- 1-sum((Gkern_SD_D_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Gkern_SD_D_UN <- 1-sum((Gkern_SD_D_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Gkern_SD_UN_D <- 1-sum((Gkern_SD_UN_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Gkern_SD_UN_UN <- 1-sum((Gkern_SD_UN_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Gkern_SD_FA_D <- 1-sum((Gkern_SD_FA_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Gkern_SD_FA_UN <- 1-sum((Gkern_SD_FA_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    
    R2b_Gkern_WI_D_D <- 1-sum((Gkern_WI_D_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Gkern_WI_D_UN <- 1-sum((Gkern_WI_D_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Gkern_WI_UN_D <- 1-sum((Gkern_WI_UN_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Gkern_WI_UN_UN <- 1-sum((Gkern_WI_UN_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Gkern_WI_FA_D <- 1-sum((Gkern_WI_FA_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Gkern_WI_FA_UN <- 1-sum((Gkern_WI_FA_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    
    #M kernel
    R2b_Mkern_MN_D_D <- 1-sum((Mkern_MN_D_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Mkern_MN_D_UN <- 1-sum((Mkern_MN_D_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Mkern_MN_UN_D <- 1-sum((Mkern_MN_UN_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Mkern_MN_UN_UN <- 1-sum((Mkern_MN_UN_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Mkern_MN_FA_D <- 1-sum((Mkern_MN_FA_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_Mkern_MN_FA_UN <- 1-sum((Mkern_MN_FA_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    
    R2b_Mkern_SD_D_D <- 1-sum((Mkern_SD_D_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Mkern_SD_D_UN <- 1-sum((Mkern_SD_D_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Mkern_SD_UN_D <- 1-sum((Mkern_SD_UN_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Mkern_SD_UN_UN <- 1-sum((Mkern_SD_UN_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Mkern_SD_FA_D <- 1-sum((Mkern_SD_FA_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_Mkern_SD_FA_UN <- 1-sum((Mkern_SD_FA_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    
    R2b_Mkern_WI_D_D <- 1-sum((Mkern_WI_D_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Mkern_WI_D_UN <- 1-sum((Mkern_WI_D_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Mkern_WI_UN_D <- 1-sum((Mkern_WI_UN_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Mkern_WI_UN_UN <- 1-sum((Mkern_WI_UN_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Mkern_WI_FA_D <- 1-sum((Mkern_WI_FA_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_Mkern_WI_FA_UN <- 1-sum((Mkern_WI_FA_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    
    #GM kernel
    R2b_GMkern_MN_D_D <- 1-sum((GMkern_MN_D_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_GMkern_MN_D_UN <- 1-sum((GMkern_MN_D_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_GMkern_MN_UN_D <- 1-sum((GMkern_MN_UN_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_GMkern_MN_UN_UN <- 1-sum((GMkern_MN_UN_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_GMkern_MN_FA_D <- 1-sum((GMkern_MN_FA_D$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    R2b_GMkern_MN_FA_UN <- 1-sum((GMkern_MN_FA_UN$yHat[tst,test.env1]-Y[tst,test.env1])^2)/sum((Y[tst,test.env1]-mean(Y[tst,test.env1]))^2)
    
    R2b_GMkern_SD_D_D <- 1-sum((GMkern_SD_D_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_GMkern_SD_D_UN <- 1-sum((GMkern_SD_D_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_GMkern_SD_UN_D <- 1-sum((GMkern_SD_UN_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_GMkern_SD_UN_UN <- 1-sum((GMkern_SD_UN_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_GMkern_SD_FA_D <- 1-sum((GMkern_SD_FA_D$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    R2b_GMkern_SD_FA_UN <- 1-sum((GMkern_SD_FA_UN$yHat[tst,test.env2]-Y[tst,test.env2])^2)/sum((Y[tst,test.env2]-mean(Y[tst,test.env2]))^2)
    
    R2b_GMkern_WI_D_D <- 1-sum((GMkern_WI_D_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_GMkern_WI_D_UN <- 1-sum((GMkern_WI_D_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_GMkern_WI_UN_D <- 1-sum((GMkern_WI_UN_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_GMkern_WI_UN_UN <- 1-sum((GMkern_WI_UN_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_GMkern_WI_FA_D <- 1-sum((GMkern_WI_FA_D$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    R2b_GMkern_WI_FA_UN <- 1-sum((GMkern_WI_FA_UN$yHat[tst,test.env3]-Y[tst,test.env3])^2)/sum((Y[tst,test.env3]-mean(Y[tst,test.env3]))^2)
    
    #--------------------------------------------------------
    #create a data frame of each run
    print("######################print predability of Gkern################")
    print(c(Predability_Gkern_MN_D_D,
            Predability_Gkern_MN_D_UN,
            Predability_Gkern_MN_UN_D,
            Predability_Gkern_MN_UN_UN,
            Predability_Gkern_MN_FA_D,
            Predability_Gkern_MN_FA_UN,
            Predability_Gkern_SD_D_D,
            Predability_Gkern_SD_D_UN,
            Predability_Gkern_SD_UN_D,
            Predability_Gkern_SD_UN_UN,
            Predability_Gkern_SD_FA_D,
            Predability_Gkern_SD_FA_UN,
            Predability_Gkern_WI_D_D,
            Predability_Gkern_WI_D_UN,
            Predability_Gkern_WI_UN_D,
            Predability_Gkern_WI_UN_UN,
            Predability_Gkern_WI_FA_D,
            Predability_Gkern_WI_FA_UN))
    print("######################print predability of Mkern################")
    print(c(Predability_Mkern_MN_D_D,
            Predability_Mkern_MN_D_UN,
            Predability_Mkern_MN_UN_D,
            Predability_Mkern_MN_UN_UN,
            Predability_Mkern_MN_FA_D,
            Predability_Mkern_MN_FA_UN,
            Predability_Mkern_SD_D_D,
            Predability_Mkern_SD_D_UN,
            Predability_Mkern_SD_UN_D,
            Predability_Mkern_SD_UN_UN,
            Predability_Mkern_SD_FA_D,
            Predability_Mkern_SD_FA_UN,
            Predability_Mkern_WI_D_D,
            Predability_Mkern_WI_D_UN,
            Predability_Mkern_WI_UN_D,
            Predability_Mkern_WI_UN_UN,
            Predability_Mkern_WI_FA_D,
            Predability_Mkern_WI_FA_UN))
    print("######################print predability of GMkern################")
    print(c(Predability_GMkern_MN_D_D,
            Predability_GMkern_MN_D_UN,
            Predability_GMkern_MN_UN_D,
            Predability_GMkern_MN_UN_UN,
            Predability_GMkern_MN_FA_D,
            Predability_GMkern_MN_FA_UN,
            Predability_GMkern_SD_D_D,
            Predability_GMkern_SD_D_UN,
            Predability_GMkern_SD_UN_D,
            Predability_GMkern_SD_UN_UN,
            Predability_GMkern_SD_FA_D,
            Predability_GMkern_SD_FA_UN,
            Predability_GMkern_WI_D_D,
            Predability_GMkern_WI_D_UN,
            Predability_GMkern_WI_UN_D,
            Predability_GMkern_WI_UN_UN,
            Predability_GMkern_WI_FA_D,
            Predability_GMkern_WI_FA_UN))
    #rbind prediction data
    results.tmp <- data.frame(Kernel=c(rep("G",18),rep("M",18),rep("G+M",18)),
                         Fold=kfold,
                         TraitName=colnames(phedat)[i],
                         Model = rep(c("D-D","D-UN","UN-D","UN-UN","FA-D","FA-UN"),9),
                         TestENV=rep(c(rep(test.env1,6), rep(test.env2,6), rep(test.env3,6)),3),
                         RunID=r,
                         PredAbility=c(Predability_Gkern_MN_D_D,
                                       Predability_Gkern_MN_D_UN,
                                       Predability_Gkern_MN_UN_D,
                                       Predability_Gkern_MN_UN_UN,
                                       Predability_Gkern_MN_FA_D,
                                       Predability_Gkern_MN_FA_UN,
                                       Predability_Gkern_SD_D_D,
                                       Predability_Gkern_SD_D_UN,
                                       Predability_Gkern_SD_UN_D,
                                       Predability_Gkern_SD_UN_UN,
                                       Predability_Gkern_SD_FA_D,
                                       Predability_Gkern_SD_FA_UN,
                                       Predability_Gkern_WI_D_D,
                                       Predability_Gkern_WI_D_UN,
                                       Predability_Gkern_WI_UN_D,
                                       Predability_Gkern_WI_UN_UN,
                                       Predability_Gkern_WI_FA_D,
                                       Predability_Gkern_WI_FA_UN,
                                       Predability_Mkern_MN_D_D,
                                       Predability_Mkern_MN_D_UN,
                                       Predability_Mkern_MN_UN_D,
                                       Predability_Mkern_MN_UN_UN,
                                       Predability_Mkern_MN_FA_D,
                                       Predability_Mkern_MN_FA_UN,
                                       Predability_Mkern_SD_D_D,
                                       Predability_Mkern_SD_D_UN,
                                       Predability_Mkern_SD_UN_D,
                                       Predability_Mkern_SD_UN_UN,
                                       Predability_Mkern_SD_FA_D,
                                       Predability_Mkern_SD_FA_UN,
                                       Predability_Mkern_WI_D_D,
                                       Predability_Mkern_WI_D_UN,
                                       Predability_Mkern_WI_UN_D,
                                       Predability_Mkern_WI_UN_UN,
                                       Predability_Mkern_WI_FA_D,
                                       Predability_Mkern_WI_FA_UN,
                                       Predability_GMkern_MN_D_D,
                                       Predability_GMkern_MN_D_UN,
                                       Predability_GMkern_MN_UN_D,
                                       Predability_GMkern_MN_UN_UN,
                                       Predability_GMkern_MN_FA_D,
                                       Predability_GMkern_MN_FA_UN,
                                       Predability_GMkern_SD_D_D,
                                       Predability_GMkern_SD_D_UN,
                                       Predability_GMkern_SD_UN_D,
                                       Predability_GMkern_SD_UN_UN,
                                       Predability_GMkern_SD_FA_D,
                                       Predability_GMkern_SD_FA_UN,
                                       Predability_GMkern_WI_D_D,
                                       Predability_GMkern_WI_D_UN,
                                       Predability_GMkern_WI_UN_D,
                                       Predability_GMkern_WI_UN_UN,
                                       Predability_GMkern_WI_FA_D,
                                       Predability_GMkern_WI_FA_UN),
                         PredAccuracy=c(Predability_Gkern_MN_D_D,
                                        Predability_Gkern_MN_D_UN,
                                        Predability_Gkern_MN_UN_D,
                                        Predability_Gkern_MN_UN_UN,
                                        Predability_Gkern_MN_FA_D,
                                        Predability_Gkern_MN_FA_UN,
                                        Predability_Gkern_SD_D_D,
                                        Predability_Gkern_SD_D_UN,
                                        Predability_Gkern_SD_UN_D,
                                        Predability_Gkern_SD_UN_UN,
                                        Predability_Gkern_SD_FA_D,
                                        Predability_Gkern_SD_FA_UN,
                                        Predability_Gkern_WI_D_D,
                                        Predability_Gkern_WI_D_UN,
                                        Predability_Gkern_WI_UN_D,
                                        Predability_Gkern_WI_UN_UN,
                                        Predability_Gkern_WI_FA_D,
                                        Predability_Gkern_WI_FA_UN,
                                        Predability_Mkern_MN_D_D,
                                        Predability_Mkern_MN_D_UN,
                                        Predability_Mkern_MN_UN_D,
                                        Predability_Mkern_MN_UN_UN,
                                        Predability_Mkern_MN_FA_D,
                                        Predability_Mkern_MN_FA_UN,
                                        Predability_Mkern_SD_D_D,
                                        Predability_Mkern_SD_D_UN,
                                        Predability_Mkern_SD_UN_D,
                                        Predability_Mkern_SD_UN_UN,
                                        Predability_Mkern_SD_FA_D,
                                        Predability_Mkern_SD_FA_UN,
                                        Predability_Mkern_WI_D_D,
                                        Predability_Mkern_WI_D_UN,
                                        Predability_Mkern_WI_UN_D,
                                        Predability_Mkern_WI_UN_UN,
                                        Predability_Mkern_WI_FA_D,
                                        Predability_Mkern_WI_FA_UN,
                                        Predability_GMkern_MN_D_D,
                                        Predability_GMkern_MN_D_UN,
                                        Predability_GMkern_MN_UN_D,
                                        Predability_GMkern_MN_UN_UN,
                                        Predability_GMkern_MN_FA_D,
                                        Predability_GMkern_MN_FA_UN,
                                        Predability_GMkern_SD_D_D,
                                        Predability_GMkern_SD_D_UN,
                                        Predability_GMkern_SD_UN_D,
                                        Predability_GMkern_SD_UN_UN,
                                        Predability_GMkern_SD_FA_D,
                                        Predability_GMkern_SD_FA_UN,
                                        Predability_GMkern_WI_D_D,
                                        Predability_GMkern_WI_D_UN,
                                        Predability_GMkern_WI_UN_D,
                                        Predability_GMkern_WI_UN_UN,
                                        Predability_GMkern_WI_FA_D,
                                        Predability_GMkern_WI_FA_UN)/sqrt(hsq),
                         R2a=c(R2a_Gkern_MN_D_D,
                               R2a_Gkern_MN_D_UN,
                               R2a_Gkern_MN_UN_D,
                               R2a_Gkern_MN_UN_UN,
                               R2a_Gkern_MN_FA_D,
                               R2a_Gkern_MN_FA_UN,
                               R2a_Gkern_SD_D_D,
                               R2a_Gkern_SD_D_UN,
                               R2a_Gkern_SD_UN_D,
                               R2a_Gkern_SD_UN_UN,
                               R2a_Gkern_SD_FA_D,
                               R2a_Gkern_SD_FA_UN,
                               R2a_Gkern_WI_D_D,
                               R2a_Gkern_WI_D_UN,
                               R2a_Gkern_WI_UN_D,
                               R2a_Gkern_WI_UN_UN,
                               R2a_Gkern_WI_FA_D,
                               R2a_Gkern_WI_FA_UN,
                               R2a_Mkern_MN_D_D,
                               R2a_Mkern_MN_D_UN,
                               R2a_Mkern_MN_UN_D,
                               R2a_Mkern_MN_UN_UN,
                               R2a_Mkern_MN_FA_D,
                               R2a_Mkern_MN_FA_UN,
                               R2a_Mkern_SD_D_D,
                               R2a_Mkern_SD_D_UN,
                               R2a_Mkern_SD_UN_D,
                               R2a_Mkern_SD_UN_UN,
                               R2a_Mkern_SD_FA_D,
                               R2a_Mkern_SD_FA_UN,
                               R2a_Mkern_WI_D_D,
                               R2a_Mkern_WI_D_UN,
                               R2a_Mkern_WI_UN_D,
                               R2a_Mkern_WI_UN_UN,
                               R2a_Mkern_WI_FA_D,
                               R2a_Mkern_WI_FA_UN,
                               R2a_GMkern_MN_D_D,
                               R2a_GMkern_MN_D_UN,
                               R2a_GMkern_MN_UN_D,
                               R2a_GMkern_MN_UN_UN,
                               R2a_GMkern_MN_FA_D,
                               R2a_GMkern_MN_FA_UN,
                               R2a_GMkern_SD_D_D,
                               R2a_GMkern_SD_D_UN,
                               R2a_GMkern_SD_UN_D,
                               R2a_GMkern_SD_UN_UN,
                               R2a_GMkern_SD_FA_D,
                               R2a_GMkern_SD_FA_UN,
                               R2a_GMkern_WI_D_D,
                               R2a_GMkern_WI_D_UN,
                               R2a_GMkern_WI_UN_D,
                               R2a_GMkern_WI_UN_UN,
                               R2a_GMkern_WI_FA_D,
                               R2a_GMkern_WI_FA_UN),
                         R2b=c(R2b_Gkern_MN_D_D,
                               R2b_Gkern_MN_D_UN,
                               R2b_Gkern_MN_UN_D,
                               R2b_Gkern_MN_UN_UN,
                               R2b_Gkern_MN_FA_D,
                               R2b_Gkern_MN_FA_UN,
                               R2b_Gkern_SD_D_D,
                               R2b_Gkern_SD_D_UN,
                               R2b_Gkern_SD_UN_D,
                               R2b_Gkern_SD_UN_UN,
                               R2b_Gkern_SD_FA_D,
                               R2b_Gkern_SD_FA_UN,
                               R2b_Gkern_WI_D_D,
                               R2b_Gkern_WI_D_UN,
                               R2b_Gkern_WI_UN_D,
                               R2b_Gkern_WI_UN_UN,
                               R2b_Gkern_WI_FA_D,
                               R2b_Gkern_WI_FA_UN,
                               R2b_Mkern_MN_D_D,
                               R2b_Mkern_MN_D_UN,
                               R2b_Mkern_MN_UN_D,
                               R2b_Mkern_MN_UN_UN,
                               R2b_Mkern_MN_FA_D,
                               R2b_Mkern_MN_FA_UN,
                               R2b_Mkern_SD_D_D,
                               R2b_Mkern_SD_D_UN,
                               R2b_Mkern_SD_UN_D,
                               R2b_Mkern_SD_UN_UN,
                               R2b_Mkern_SD_FA_D,
                               R2b_Mkern_SD_FA_UN,
                               R2b_Mkern_WI_D_D,
                               R2b_Mkern_WI_D_UN,
                               R2b_Mkern_WI_UN_D,
                               R2b_Mkern_WI_UN_UN,
                               R2b_Mkern_WI_FA_D,
                               R2b_Mkern_WI_FA_UN,
                               R2b_GMkern_MN_D_D,
                               R2b_GMkern_MN_D_UN,
                               R2b_GMkern_MN_UN_D,
                               R2b_GMkern_MN_UN_UN,
                               R2b_GMkern_MN_FA_D,
                               R2b_GMkern_MN_FA_UN,
                               R2b_GMkern_SD_D_D,
                               R2b_GMkern_SD_D_UN,
                               R2b_GMkern_SD_UN_D,
                               R2b_GMkern_SD_UN_UN,
                               R2b_GMkern_SD_FA_D,
                               R2b_GMkern_SD_FA_UN,
                               R2b_GMkern_WI_D_D,
                               R2b_GMkern_WI_D_UN,
                               R2b_GMkern_WI_UN_D,
                               R2b_GMkern_WI_UN_UN,
                               R2b_GMkern_WI_FA_D,
                               R2b_GMkern_WI_FA_UN),
                         stringsAsFactors = F)
    #output results of a single run for each trait
    write.csv(results.tmp,
              paste0(DIR_output_SingleRun,
                     "ElitePanel_multiomics_prediction_runID-",r,"_",colnames(phedat)[i],".csv"))
    
    return(results.tmp)
    results.tmp <- NULL #empty results.tmp
    
  }#end of foreach, 2nd layer loop
  
  #output all results of trait i
  write.csv(results,
            paste0(DIR_output,"ElitePanel_multiomics_prediction_TraitID-",i,"_",colnames(phedat)[i],"_",kfold,"fold_","rundate",rundate,".csv"),
            row.names = F)
  
  print(paste("PhenoID =",i,colnames(phedat)[i],"is done!"))
  results = NULL #empty
  
}#end of first layer loop

#shut down clusters
stopCluster(cl)
proc.time() - ptm # Stop the clock
