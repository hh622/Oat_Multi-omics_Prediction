#=================================================
#script01b_DivPanel_prepare_relationship_matrices
#HH 2020-07-30
#=================================================

#----------------------------------------
# 2021-07-30
# prepare_relationship_matrices


#clear workspace
rm(list=ls()) 

#loading libraries
library(rrBLUP)
library(sommer)
library(lineup)
library(tibble) #to load add_column()
library(lineup)
# source('../../BSFG/Estimate_gcor_prediction.R')

#--------------------------------
#SECTION 1: loading omics data
#--------------------------------
#set dir
DIR_omicsdat="../01_Omicsdat/"
list.files(DIR_omicsdat)

#loading omics data
ls1=load(file = paste0(DIR_omicsdat,"DivPanel_Omicsdata.Rdata"))
ls1# [1] "phedat"   "genodat"  "transdat" "gcmsdat"  "lcmsdat" 

#checking whether all omics data have the same set of genotypes
sapply(list(rownames(genodat), rownames(transdat), rownames(gcmsdat), rownames(lcmsdat)), 
       FUN = identical, rownames(phedat)) #all returned TRUE

#checking dimensiosn of all omics data
sapply(list(phedat, genodat, transdat, gcmsdat, lcmsdat), 
       FUN = dim) #all returned TRUE
#      [,1]  [,2]  [,3] [,4] [,5]
# [1,]  333   333   333  333  333
# [2,]   19 73014 54169  793 1233

#loading heritability of omics data
ls2=load(file=paste0(DIR_omicsdat,"DivPanel_Omicsdata_heritability.Rdata"))
ls2 #[1] "hsq_Div_pheno"  "hsq_Div_rnaseq" "hsq_Div_gcms"   "hsq_Div_lcms"

#-------------------------------------------------
#SECTION 2: calculating relationship matrices
#-------------------------------------------------

#--------------------------------------------------
#step 2-1: calculate genomic relationship matrix(GRM)
genodat[1:3,1:4]# genotypes in row, markers in column
# checking dim
dim(genodat) #333 73014 
# checking coding
unique(c(as.matrix(genodat)))# -1  1  0

#prep genotype matrix
genodat <- genodat[,apply(genodat,2,var)!=0] 
genodat <- round(genodat,0)

# Calculates the realized additive relationship matrix
Amat_rrblup <- rrBLUP::A.mat(as.matrix(genodat))
Amat_sommer <- sommer::A.mat(as.matrix(genodat))
all.equal(Amat_rrblup, Amat_sommer) #TRUE
Amat <- Amat_rrblup

# Calculates the realized dominance relationship matrix
Dmat <- sommer::D.mat(as.matrix(genodat))

# Calculates the realized epistatic relationship matrix of second order
Emat <- sommer::E.mat(as.matrix(genodat))

#checking
sapply(list(Amat, Dmat, Emat), FUN = dim) 
#      [,1] [,2] [,3]
# [1,]  333  333  333
# [2,]  333  333  333
sapply(list(rownames(Amat), rownames(Dmat), rownames(Emat), 
            colnames(Amat), colnames(Dmat), colnames(Emat)), 
       FUN = identical, rownames(genodat))#all returned TRUE

#-----------------------------------------------------------
#step 2-2: calculate transcriptomic relationship matrix(TRM)
dim(transdat)#333 54169
sum(is.na(transdat)) #24642

trans1 <- transdat[,hsq_Div_rnaseq$hsq>0 & !is.na(hsq_Div_rnaseq$hsq)];   dim(trans1) #333 38766
trans2 <- transdat[,hsq_Div_rnaseq$hsq>0.2 & !is.na(hsq_Div_rnaseq$hsq)]; dim(trans2) #333 26058
trans3 <- transdat[,hsq_Div_rnaseq$hsq>0.4 & !is.na(hsq_Div_rnaseq$hsq)]; dim(trans3) #333 12637
TRM1 <- (as.matrix(trans1) %*% as.matrix(t(trans1)))/ncol(trans1)
TRM2 <- (as.matrix(trans2) %*% as.matrix(t(trans2)))/ncol(trans2)
TRM3 <- (as.matrix(trans3) %*% as.matrix(t(trans3)))/ncol(trans3)

hist(corbetw2mat(TRM1, TRM2))
hist(corbetw2mat(TRM1, TRM3))
#Notes: there is high correlation using hsq>0, hsq >0.2, hsq>0.4 (r>=.99)

#decide TRM
TRM <- TRM3

#checking
dim(TRM)#333 333
TRM[1:3,1:4]

#----------------------------------------------------------
#step 2-3: calculate metabolomic relationship matrix (MRM)

#merge gcms + lcms
if(identical(rownames(gcmsdat),rownames(lcmsdat))){
  metadat <- cbind(gcmsdat, lcmsdat)
}
dim(metadat)#333 2026

#merge hsq of gcms + lcms
if(identical(colnames(hsq_Div_gcms),colnames(hsq_Div_lcms))){
  hsq_Div_meta <- rbind(hsq_Div_gcms, hsq_Div_lcms)
}
head(hsq_Div_meta)

#truncate metadata according to hsq
if(identical(colnames(metadat),hsq_Div_meta$Traitname))
{
  metadat1 <- metadat[,hsq_Div_meta$hsq>0]
  metadat2 <- metadat[,hsq_Div_meta$hsq>0.2]
  metadat3 <- metadat[,hsq_Div_meta$hsq>0.4]
}
sapply(list(metadat, metadat1, metadat2, metadat3), FUN = dim) 
#     [,1] [,2] [,3] [,4]
# [1,]  333  333  333  333
# [2,] 2026 1696 1377  861

#Calculate MRMs
MRM1 <- (as.matrix(metadat1) %*% as.matrix(t(metadat1)))/ncol(metadat1)
MRM2 <- (as.matrix(metadat2) %*% as.matrix(t(metadat2)))/ncol(metadat2)
MRM3 <- (as.matrix(metadat3) %*% as.matrix(t(metadat3)))/ncol(metadat3)
hist(corbetw2mat(MRM1, MRM2))
hist(corbetw2mat(MRM1, MRM3))
#Notes: there is high correlation using hsq>0, hsq >0.2, hsq>0.4 (r>=.95)

#decide MRM
MRM <- MRM1

#---------------------------------
#step 2-4: save data
sapply(list(rownames(Amat),colnames(Amat),
            rownames(Dmat),colnames(Dmat),
            rownames(Emat),colnames(Emat),
            rownames(TRM),colnames(TRM),
            rownames(MRM),colnames(MRM),
            rownames(genodat)), 
       FUN = identical, 
       rownames(phedat)) #all returned TRUE
sapply(list(Amat, Dmat, Emat, TRM, MRM), FUN = dim) 
#      [,1] [,2] [,3] [,4] [,5]
# [1,]  333  333  333  333  333
# [2,]  333  333  333  333  333

sapply(list(metadat, metadat1, metadat2, metadat3), FUN = dim) 

#set dir
DIR_output <- "../01_Omicsdat/"
if(!dir.exists(DIR_output)) dir.create(DIR_output)
list.files(DIR_output)

#add phedat and trait hsq to the multiomics prediction input data
head(hsq_Div_pheno)
head(phedat)
identical(colnames(phedat), hsq_Div_pheno$Traitname)

#save relationship matrix data
save(Amat, Dmat, Emat, TRM, MRM,
     file = paste0(DIR_output,"DivPanel_Omicsdata_GRM_TRM_MRM_333x333.RData"))

# save pheno data
save(phedat, hsq_Div_pheno,
     file = paste0(DIR_output,"DivPanel_pheno.RData"))




