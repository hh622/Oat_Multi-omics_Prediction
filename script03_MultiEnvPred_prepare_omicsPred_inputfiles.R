#==================================================
#script01_ElitePanel_prepare_omicsPred_input
#HH 2021-08-01
#==================================================

#notes 2021-08-01
#after discussion with Dan (2020-08-20), HH decided to use the other two-ENV to make BULPs to 
#predict the target phenotypes of the 3rd ENV

#clear workspace
rm(list=ls()) 

# load libraries
library(rrBLUP)
library(lineup) #corbetw2mat

#set output dir
DIR_Output <- "../01_Omicsdat/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

#-------------------------
#SECTION 01: read in data
#-------------------------

# DIR_Pheno <- "/workdir/hh622/P2/25_Elite17_Phenotypes_all/2_output/"
# DIR_Geno <- "/workdir/hh622/P2/03_genodat_RNAseq18/"
# DIR_meta <- "/workdir/hh622/P2/25_Elite17_Phenotypes_all/2_output/"

DIR_Pheno <- "/home/haixiao/workdir/P2/25_Elite17_Phenotypes_all/2_output/"
DIR_Geno <- "/home/haixiao/workdir/P2/03_RNAseq18_Genotypes/"
DIR_GC02LC02_2ENV <- "/home/haixiao/workdir/P2/25_Elite17_Phenotypes_all/2_output/2ENVs_hsq_BLUPs_GC02LC02/"

list.files(DIR_Pheno)
list.files(DIR_Geno)
list.files(DIR_GC02LC02_2ENV)

#--------------------------------------
#SECTION 02: Prepare Phenotypic data
#--------------------------------------

#-------------------------
#read in phenotypic data
phedat0 <- read.csv(paste0(DIR_Pheno,"Elite17_Phenodat_LC02_BLUPs_within_site.csv"),
                    header=T, stringsAsFactors = F)

head(phedat0)
dim(phedat0) #673  28

#----------------------------------------------------
#read in broad-sense heritability of RNASeq18 traits
trait.H2.singlesite <- read.csv(paste0(DIR_Pheno,"Elite17_Phenodat_LC02_hsq_within_site.csv"),
                                header = T,stringsAsFactors = F)
head(trait.H2.singlesite)
dim(trait.H2.singlesite)#60  3

trait.H2.singlesite[trait.H2.singlesite$hsq < 0.2,]
# LOC           Traitname       hsq
# 23  SD          Seedlength 0.1972382
# 24  SD           Seedwidth 0.0000000
# 47  WI hundred_hull_weight 0.1415897

trait.H2.acrosssite <- read.csv(paste0(DIR_Pheno,"Elite17_Phenodat_LC02_hsq_across_sites.csv"),
                                header = T,stringsAsFactors = F)

#------------------------------------------------------
#exclude three traits that have low hsq in single site
phedat0 <- phedat0[,!(colnames(phedat0)%in%(trait.H2.singlesite$Traitname[trait.H2.singlesite$hsq < 0.2]))]
dim(phedat0)#[1] 673  25
head(phedat0)
trait.H2.acrosssite <- trait.H2.acrosssite[!(trait.H2.acrosssite$Traitname %in% (trait.H2.singlesite$Traitname[trait.H2.singlesite$hsq < 0.2])),]
trait.H2.singlesite <- trait.H2.singlesite[!(trait.H2.singlesite$Traitname%in%(trait.H2.singlesite$Traitname[trait.H2.singlesite$hsq < 0.2])),]

dim(trait.H2.singlesite)
dim(trait.H2.acrosssite)
head(trait.H2.singlesite)
head(trait.H2.acrosssite)

#-------------------------------------------
#SECTION 03: Prepare genotypic data 
#-------------------------------------------

#-----------------------
#read in genotypic data
list.files(DIR_Geno)
genodat0 <- read.table(paste0(DIR_Geno,"SNPs_genomat_72719x568_glmimputed.txt"),header = T,stringsAsFactors = F)

dim(genodat0)
sum(rownames(genodat0) %in% phedat0$LINE)#230
genodat0[1:3,1:4]
setdiff(phedat0$LINE,rownames(genodat0))
#"IL10-9867" "IL00-654"

#-------------------------------------------
#SECTION 04: Prepare metabolomic data
#-------------------------------------------

# Load the data saved in the first part
ls=load(file = paste0(DIR_GC02LC02_2ENV,"Elite17_GC02qc_plus_LC02qc.RData"))
ls
# [1] "blup.gc02lc02.Elite17qc_MN_SD"  "blup.gc02lc02.Elite17qc_MN_WI" 
# [3] "blup.gc02lc02.Elite17qc_SD_WI"  "hsq.gc02lc02.Elite17qc_MN_SD2b"
# [5] "hsq.gc02lc02.Elite17qc_MN_WI2b" "hsq.gc02lc02.Elite17qc_SD_WI2b"

mat1=as.matrix(blup.gc02lc02.Elite17qc_MN_SD)
mat2=as.matrix(blup.gc02lc02.Elite17qc_MN_WI)
mat3=as.matrix(blup.gc02lc02.Elite17qc_SD_WI)
sum(is.na(mat1)); sum(is.na(mat2)); sum(is.na(mat3))#0; 0; 0
#There is no NA left

# ##use median to substitute NAs
# metadat0_MN_SD <- sweep(mat1, MARGIN = 2, 
#               STATS = apply(mat1, 2, median, na.rm=TRUE),
#               FUN =  function(x,s) ifelse(is.na(x), s, x)
# )
# metadat0_MN_WI <- sweep(mat2, MARGIN = 2, 
#               STATS = apply(mat2, 2, median, na.rm=TRUE),
#               FUN =  function(x,s) ifelse(is.na(x), s, x)
# )
# metadat0_SD_WI <- sweep(mat3, MARGIN = 2, 
#               STATS = apply(mat3, 2, median, na.rm=TRUE),
#               FUN =  function(x,s) ifelse(is.na(x), s, x)
# )
metadat0_MN_SD <- mat1
metadat0_MN_WI <- mat2
metadat0_SD_WI <- mat3

metadat0_MN_SD[1:3,1:4]
sum(is.na(metadat0_MN_SD))#0
sum(is.na(metadat0_MN_WI))#0
sum(is.na(metadat0_SD_WI))#0
str(metadat0_MN_SD)

#------------------------------------------------------------------
#SECTION 05: truncate datasets to make them all have them LINEs
#------------------------------------------------------------------

#extract those lines with pheotypic values in all three ENVs
phedat <- phedat0[phedat0$LINE %in% (names(table(phedat0$LINE))[table(phedat0$LINE)==3]),]

#get a set of common lines
lines.com <- 
  Reduce(intersect, list(phedat$LINE,
                         rownames(genodat0),
                         rownames(metadat0_MN_SD),
                         rownames(metadat0_MN_WI),
                         rownames(metadat0_SD_WI)))
length(lines.com)

#truncate datasets
phedat <- phedat[phedat$LINE %in% lines.com,]
genodat <- genodat0[rownames(genodat0) %in% lines.com,]
metadat_MN_SD <- metadat0_MN_SD[rownames(metadat0_MN_SD) %in% lines.com,]
metadat_MN_WI <- metadat0_MN_WI[rownames(metadat0_MN_WI) %in% lines.com,]
metadat_SD_WI <- metadat0_SD_WI[rownames(metadat0_SD_WI) %in% lines.com,]

dim(phedat)#651  25
dim(genodat) #217 72719
dim(metadat_MN_SD)#217 890
dim(metadat_MN_WI)#217 806
dim(metadat_SD_WI)#217 778
#=> so the truncated data has 217 genotypes

#make genotypic and metabolomic datasets have the same order in LINE
phedat <- phedat[,!(colnames(phedat) %in% c("BLK","ROW" ,"COL","PLOT","BATCH.FAME","BATCH.LC02"))]
genodat <- genodat[rownames(metadat_MN_SD),] #re-order rows of genodat

identical(setdiff(phedat$LINE,rownames(genodat)),character(0)) #TRUE
identical(setdiff(rownames(genodat),phedat$LINE),character(0)) #TRUE

sapply(list(rownames(metadat_MN_SD),
            rownames(metadat_MN_WI),
            rownames(metadat_SD_WI)),
       FUN = identical,
       rownames(genodat)) #TRUE TRUE TRUE
# => so LINEs in genotypic and metabolomic datasets have the same order

#----------------------------------------
#SECTION 06: calculate GRM and MRM
#----------------------------------------

#------------------------
#setp 6-1: calculate GRM

#round imputed genotype value to a nearest integer
genodat <- round(genodat)

# calculate genomic relationship matrix
dim(genodat)# 217 72719
Ksnp <- A.mat(genodat) 
Ksnp[1:3,1:4]
dim(Ksnp)#217 217

#--------------------------
#step 6-2: calculate MRM
#e.g. metadat <- as.matrix(metadat);Kmeta <- (metadat%*% t(metadat))/ncol(metadat)

KM_MN_SD <- (metadat_MN_SD %*% t(metadat_MN_SD))/ncol(metadat_MN_SD)
KM_MN_WI <- (metadat_MN_WI %*% t(metadat_MN_WI))/ncol(metadat_MN_WI)
KM_SD_WI <- (metadat_SD_WI %*% t(metadat_SD_WI))/ncol(metadat_SD_WI)

dim(KM_MN_SD)#217 217
dim(KM_MN_WI)#217 217
dim(KM_SD_WI)#217 217
KM_MN_SD[1:3,1:4]
KM_MN_WI[1:3,1:4]
KM_SD_WI[1:3,1:4]


#----------------------------------------
#step 6-3: output GRM and MRM
sapply(list(rownames(Ksnp),colnames(Ksnp),
            rownames(KM_MN_SD),colnames(KM_MN_SD),
            rownames(KM_MN_WI),colnames(KM_MN_WI),
            rownames(KM_SD_WI),colnames(KM_SD_WI)), 
       FUN = identical, 
       rownames(Ksnp)) #all returned TRUE
#notes: for phedat at each LOC, the order of LINE is not yet identical to GRM/MRM
sum(is.na(Ksnp)) #0
sum(is.na(KM_MN_SD)) #0
sum(is.na(KM_MN_WI)) #0
sum(is.na(KM_SD_WI)) #0

identical(colnames(phedat)[-1*1:2],trait.H2.acrosssite$Traitname)

save(phedat,
     trait.H2.acrosssite,
     trait.H2.singlesite,
     Ksnp,
     KM_MN_SD, KM_MN_WI, KM_SD_WI,
     file = paste0(DIR_Output,"ElitePanel_OmicsPred_Input_GC02LC02.RData"))

ls2=load(file = paste0(DIR_Output,"ElitePanel_OmicsPred_Input_GC02LC02.RData"))
ls2
