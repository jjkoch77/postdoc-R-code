##############################################################################################
######## PART I: OXBS-MLE PREP AND PROCESSING ###########
##############################################################################################

## 1. Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("minfi", version = "3.8")
BiocManager::install("IlluminaHumanMethylation450kmanifest")
BiocManager::install("IlluminaHumanMethylation450kanno.ilmn12.hg19")
BiocManager::install("ChAMP", version = "3.8")
BiocManager::install("ENmix", version = "3.8")
BiocManager::install("tidyr", version = "3.8")
BiocManager::install("dplyr", version = "3.8")
BiocManager::install("GEOquery")

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(ChAMP)
library(ENmix)
library(tidyr)
library(dplyr)
library(GEOquery)

## 2. Prep sample sheets

#Using publicly available AD BS/oxBS-450K data from the following GEO database entry:
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105109

#Phenotype data for the above study is available as SOFT files:
gds <- getGEO("GSE105109") #Large list with 7 slots
gds.pheno <- gds$GSE105109_series_matrix.txt.gz@phenoData #Pulls out annotated data frame of phenotype data
gds.pheno.df <- gds.pheno@data #Pulls out annotated data frame of phenotype data

#Remove the cerebellum samples from the phenotype data (rows 193-384):
gds.pheno.df1 <- gds.pheno.df[-c(193:384), ]

#Write.csv of the phenotype data:
write.csv(gds.pheno.df1, "/path/to/GSE105109.pheno.052819.csv")

# If phenotypic/meta data has already been prepped in a previous session, it can be directly imported here.
samples_data <- read.csv( "/path/to/GSE105109.pheno.052819.csv", header = TRUE) 

#Now try subsetting the samples data by several variables based on study design.
#First split by sex:
samples_data.m <- samples_data[(samples_data$gender.ch1 == "M"), ] #108 observations
samples_data.f <- samples_data[(samples_data$gender.ch1 == "F"), ] #84 observations

#Now subset male and female data down to control group and braak stage VI
samples_data.m.1 <- samples_data.m[(samples_data.m$post.mortem.diagnosis.ch1 == "Control") | (samples_data.m$braak.stage.ch1 == "VI"), ] #72 observations; n=36 -- n = 14 control, n = 22 AD
samples_data.f.1 <- samples_data.f[(samples_data.f$post.mortem.diagnosis.ch1 == "Control") | (samples_data.f$braak.stage.ch1 == "VI"), ] #60 observations; n=30 -- n = 14 control, n = 16 AD

#Due to slightly larger sample size, for downstream analysis, use the males (Control vs. Braak stage VI).
samples_data.m.1 <- samples_data.m[(samples_data.m$post.mortem.diagnosis.ch1 == "Control") | (samples_data.m$braak.stage.ch1 == "VI"), ] #72 observations; n=36 -- n = 14 control, n = 22 AD

#write.csv(samples_data.m.1, "/path/to/GSE105109.pheno.male.filtered.052819.csv")

#Assuming the phenotype data has already been processed before, read it in:
samples_data.m.1 <- read.csv("/path/to/GSE105109.pheno.male.filtered.052819.csv")
#Note: The version saved on the local drive has added columns for Basename, Array, and Assay (added manually in excel)

## 3. Import idat files and set annotation
#Use minfi to import idat files separated by tissue type

# Use minfi (read.metharry.exp) to import idat files using the targets defined in phenotype data frame.
# The class RGSet is a RGChannelSet object. This is the initial object of a 'library(minfi)' analysis that contains the raw intensities in the green and red channels. 
# Keep only the probes that are consistent across all arrays
# Illumina updated the EPIC recently and tossed out a couple thousand probes
# Thus, force=TRUE will only keep probes found on all arrays

rgset.m <- read.metharray.exp(targets = samples_data.m.1, recursive = TRUE, verbose = TRUE, force = TRUE)

#Can also view the manifest
manifest <- getManifest(rgset.m)
manifest

rgset.m@annotation = c(array='IlluminaHumanMethylation450k', annotation='ilmn12.hg19')

# Save rgsets for easy loading if required
save(rgset.m, file = "/path/to/GSE105109_rgsets.RData")

# If rgsets have already been generated in a previous session, they can be directly imported here.
load("/path/to/GSE105109_rgsets.RData")


## 4. Check control probes
#Use plotCtrl in the ENmix pacakge to assess internal control probes to evaluate if the chip worked properly.

#Internal Control Probes Graphs
dir.create("/path/to/Control_Probes/")
setwd("/path/to/Control_Probes/")
#rgset.all <- read.metharray.exp(targets = meta, recursive = TRUE, verbose = FALSE, force = TRUE)
dim(rgset.m)
rgset.m@annotation = c(array='IlluminaHumanMethylation450k', annotation='ilmn12.hg19')
plotCtrl(rgset.m) 
#### Plan to include these images in PDF output####

## 5. Generate MethylSet
# Use minfi to generate a MethylSet containing methylated and unmethylation signals with preprocessRaw function. 
# This calculates values without normalization.

# Set sample names
sampleNames(rgset.m) = rgset.m[[2]]

#Pull out phenotype data to examine
pd.m <- pData(rgset.m)
pd.m[,1:6]

# Set slide column to factor to avoid issues with champ.SVD
pd.m$Slide <- as.factor(pd.m$Slide) 

# Generate a MethylSet
mset.m = preprocessRaw(rgset.m)
#MSet has 485512 elements

# Generate detection p-values at the probe level across samples
detP.m <- detectionP(rgset.m)

# Extract the beta values
raw_betas.m <- getBeta(mset.m, "Illumina")

# Save msets
# Save msets for easy loading if required
save(mset.m, file = "/path/to/GSE105109_msets.RData")
save(pd.m, detP.m, raw_betas.m, file = "/path/to/GSE105109_other_raw_data.Rdata")

# Import msets
# If msets have already been generated in a previous session, they can be directly imported here.
load("/path/to/GSE105109_msets.RData")
load("/path/to/GSE105109_other_raw_data.Rdata")


## 6. Filter bad probes and/or bad samples by detection p-value
# Filter out low quality probes based on detection value
raw_betas.m[detP.m >= 0.01] <- NA

# Calculate the proportion of probes that failed the detection p-value threshold
numfail.m <- matrix(colMeans(is.na(raw_betas.m)))
numfail.m #Small percentage of probes failed across samples.

# Rename rows and columns
rownames(numfail.m) <- colnames(detP.m)
colnames(numfail.m) <- "Failed CpG Fraction"
print(numfail.m)

# Identify how many samples have less than 0.1 fraction of total probes
RemainSample.m <- which(numfail.m < 0.1)
RemainSample.m
# One sample ID -- sample 50 (GSM2818137) fails by this standard. Remove this sample and its paired data from the list.
RemainSample.m.2 <- RemainSample.m[-49] #Remove paired data -- sample 49 -- for sample 50 (failed sample)
RemainSample.m.2


#Filter the data keeping only the samples that pass with a reasonable fraction of failing probes
rgset.m <- rgset.m[, RemainSample.m.2]
detP.m <- detP.m[, RemainSample.m.2]
mset.m <- mset.m[, RemainSample.m.2]
pd.m <- pd.m[RemainSample.m.2, ]
raw_betas.m <- raw_betas.m[, RemainSample.m.2]

#Filter the remaining samples using p-value cutoff. Failed probes have a detection p-value less than 0.05.
#Set the probe cutoff to drop probes that failed in ANY sample (i.e. they must be present in all samples).
#NOTE: May be more appropriate to adjust these values (e.g. p-value > 0.01 in >10% of samples); this will vary depending on sample size.

dim(mset.m)
# [1] 485512     70

ProbeCutoff <- 0.05
mset.m.f <- mset.m[rowSums(is.na(raw_betas.m)) <= ProbeCutoff * ncol(detP.m), ]
raw_betas.m.f <- raw_betas.m[rowSums(is.na(raw_betas.m)) <= ProbeCutoff * ncol(detP.m), ]
dim(mset.m.f)
# [1] 464137     70
485512-464137
# 21375 probes removed

# Save filtered msets and raw_betas for easy loading if required
save(mset.m.f, file = "/path/to/GSE105109_msets_f.RData")
save(raw_betas.m.f, file = "/path/to/GSE105109_raw_betas_f.Rdata")

# If msets and raw betas have already been filtered in a previous session, they can be directly imported here.
load("./RData/msets_f.RData")
load("./RData/raw_betas_f.Rdata")


## 7. Mask probes from Chen et al. 2013's lists of cross-reactive probes, common SNPs at CpG probes, and SBE SNPs for hg19

###List of cross-reactive/non-specific probes to be removed
cross.react <- read.csv('/path/to/48639-non-specific-probes-Illumina450k.csv', head = T, as.is = T)
cross.react.probes <- as.character(cross.react$TargetID)

###List of probes with SNPs at CpG positions that will be excluded. 
commonSNP = read.table("/path/to/48640-polymorphic-CpGs-Illumina450k_CpG.txt", header = TRUE,fill=TRUE)##This is a list with all SNPs anywhere in the probe and with a wide variety of MAF
##Fine-tune the common SNPs list (only those with MAF>0.05):
commonSNP.f <- commonSNP[commonSNP$AF > 0.05, ]

###List of probes with SNPs at SBE that will be excluded. 
SNP.SBE = read.table("/path/to/48640-polymorphic-CpGs-Illumina450k_SBE.txt", header = TRUE,fill=TRUE)##This is a list with all SNPs anywhere in the probe and with a wide variety of MAF
##Fine-tune the common SNPs list (only those with MAF>0.05):
SNP.SBE.f <- SNP.SBE[SNP.SBE$AF > 0.05, ]


NSP_SNP = c(as.vector(cross.react$TargetID),as.vector(commonSNP.f$PROBE),as.vector(SNP.SBE.f$PROBE))##Combines both lists of cross-reactive and SNPs
NSP_SNP_nodup = unique(NSP_SNP) ##removes duplicates so we can see how many unique probes will be removed -- 81556 probes

#Now remove the cross-reactive and SNP probes from the data:
mset.m.fm <- mset.m.f[!featureNames(mset.m.f) %in% NSP_SNP_nodup, ]
raw_betas.m.fm <- raw_betas.m.f[!rownames(raw_betas.m.f) %in% NSP_SNP_nodup, ]
dim(mset.m.fm)
#[1] 386245     70
# 81556 probes removed

# Save masked msets and raw_betas for easy loading if required
save(mset.m.fm, file = "/path/to/GSE105109_msets_fm.RData")
save(raw_betas.m.fm, file = "/path/to/GSE105109_raw_betas_fm.Rdata")

# If msets and raw betas have already been masked in a previous session, they can be directly imported here.
load("./RData/msets_fm.RData")
load("./RData/raw_betas_fm.Rdata")



## 8. Adjust beta-values
#Fix Beta values that are either 0 or greater than or equal to 1
#THIS IS REALLY IMPORTANT FOR DOWNSTREAM TRANSFORMATIONS
if (min(raw_betas.m.fm, na.rm = TRUE) <= 0)
  raw_betas.m.fm[raw_betas.m.fm <= 0] <- min(raw_betas.m.fm[raw_betas.m.fm > 0])
message("Zeros in your dataset have been replaced with smallest positive value.\n")

if (max(raw_betas.m.fm, na.rm = TRUE) >= 1)
  raw_betas.m.fm[raw_betas.m.fm >= 1] <- max(raw_betas.m.fm[raw_betas.m.fm < 1])
message("Ones in your dataset have been replaced with largest value below 1.\n")


## 9. Extract raw data and run QC
#Get the intensity values and detection p-values using minfi to feed to QC in ChAMP
intensity.m <- minfi::getMeth(mset.m.fm) + minfi::getUnmeth(mset.m.fm)
detP.m <- detP.m[which(row.names(detP.m) %in% row.names(raw_betas.m.fm)), ]

#Compile the data into a list object to feed into ChAMP QC
preprocessed.raw.data.m <- list(mset = mset.m.fm, rgSet = rgset.m, pd = pd.m, intensity = intensity.m, beta = raw_betas.m.fm, detP = detP.m)

#Run a quick QC with ChAMP
champ.QC(beta = preprocessed.raw.data.m$beta, pheno = pd.m$post.mortem.diagnosis.ch1, resultsDir = "/path/to/GSE105109_CHAMP_QC")
#Obvious split between 5-mC and 5-hmC data. No obvious clustering by control vs. alzheimers within those two categories.


## 10. Dye bias correction with ssNoob (within minfi) and run QC
#Use ssNoob in minfi to normalize data

# Perform dye bias correction for beta values with ssNoob on probes remaining after filtering and masking. We can leave BS and oxBS together because this is single sample Noob.
mset.m.n <- preprocessNoob(rgset.m, dyeMethod = "single")[rownames(raw_betas.m.fm), ]
head(mset.m.n)
mset.m.n
mset.m.n[[1]]
mset.m.n[[2]]
mset.m.n[[3]]


# Set sample names
sampleNames(mset.m.n) = mset.m.n[[2]]

# Extract data from mset
betas.m.n <- getBeta(mset.m.n, "Illumina")
betas.m.n

# Compile the data into a list object to feed into ChAMP QC but with dye bias corrected beta values to re-run QC
ssNoob.data.m <- list(mset = mset.m.n, rgSet = rgset.m, pd = pd.m, intensity = intensity.m, beta = betas.m.n, detP = detP.m)

champ.QC(beta = ssNoob.data.m$beta, pheno = pd.m$post.mortem.diagnosis.ch1, resultsDir = "/path/to/GSE105109_CHAMP_QC_dyebiascorr")
champ.QC(beta = ssNoob.data.m$beta, pheno = pd.m$Assay, resultsDir = "/path/to/GSE105109_CHAMP_QC_dyebiascorr_2")
#Obvious splitting by 5-mC vs. 5-hmC, but no obvious clustering by control vs. alzheimers within those two categories.

# Perform SVD to check for batch effects
champ.SVD(beta = betas.m.n, pd = pd.m, resultsDir = "/path/to/GSE105109_CHAMP_SVD/")
#Obvious significant effect of assay (5-mC vs. 5-hmC), which is unsurprising. 
#Also, some other potential covariates to consider including in the model!

# Save normalized values}
# Save normalized values without duplicates for easy loading if required
save(betas.m.n, mset.m.n, file = "/path/to/GSE105109_ssNoob_data.RData")

# If normalized values have already been generated in a previous session, they can be directly imported here.
load("/path/to/GSE105109_ssNoob_data.RData")


## 11. Estimate cell type proportions with CETS
# Use betas corrected by ssNoob with bad probes filtered out, probes maskes, only BS (1) assay
# betas.m.n

#Remove the two samples that were filtered out for quality reasons from the phenotype ID meta data sheet, then read in new version:
samples_data.m.2 <- read.csv("/path/to/GSE105109.pheno.male.filtered.v2.052819.csv")


# Select BS Asssay (1) from these
bs.m <- as.character(samples_data.m.2$X[which(samples_data.m.2$Assay == 1)])
betas.m.n.a <- betas.m.n[,bs.m]

#Load CETS
load("/path/to/CETS_3.03/CETS_Image.RData")

# Use glial cell proportion, control brains, Caucasian samples and male samples to generate reference 
idx <- list(controlNeuron = pdBrain$celltype == "N" & pdBrain$diag == "Control" & pdBrain$ethnicity == "Caucasian" & pdBrain$sex == "Male", controlGlia = pdBrain$celltype == "G" & pdBrain$diag == "Control" & pdBrain$ethnicity == "Caucasian" & pdBrain$sex == "Male")
refProfile <- getReference(brain, idx)
head(refProfile)

# Estimate neuronal proportion
prop.m <- estProportion(betas.m.n.a, profile = refProfile)
round(prop.m, 3)

# Convert to glial proportion (1- neuronal proportion)
prop.m.g <- as.data.frame(1 - prop.m) 

# Change rowname to column
prop.m.g <- tibble::rownames_to_column(prop.m.g, var = "rowname")

# Change column name
names(prop.m.g)[1] <- "Sample"

# change column name to glial
names(prop.m.g)[2] <- "glial"

# Add cell proportion data to meta data
samples_data.m.2$glial <- NA #creates empty column for glial variable
samples_data.m.2$glial <- prop.m.g$glial[match(samples_data.m.2$X, prop.m.g$Sample)]

# Add cell proportion data to pd.c and pd.p
pd.m@listData$glial <- NA
pd.m@listData$glial <- prop.m.g$glial[match(pd.m@listData$X, prop.m.g$Sample)]

# Save new meta data as text file
write.csv(samples_data.m.2, "/path/to/GSE105109.pheno.male.filtered.v3.052919.csv")

# If meta data has already been prepped in a previous session, it can be directly imported here.
samples_data.m.3 <- read.csv("/path/to/GSE105109.pheno.male.filtered.v3.052919.csv", header = TRUE)

# Save new phenotypic data with glial proportions d
save(pd.m, file = "/path/to/GSE105109_pd_m.RData")

#Load updated phenotypic data
load("/path/to/GSE105109_pd_m.RData")


## 13. Generate full beta value data frame prior to OxBS MLE
# Create beta value and intensity matrices required for further OXBS MLE and differential methylation/hydroxymethylation analysis
samples_data.m.3 <- read.csv("/path/to/GSE105109.pheno.male.filtered.v3.052919.csv", header = TRUE)
full.m <- as.character(samples_data.m.3[which(samples_data.m.3$Replicate == 1), "X"]) #Pull out sample names ("X")
betas.m.n.full <- betas.m.n[,full.m]
intensity.m.full <- intensity.m[,full.m]


#Save data with duplicates removed
# Save normalized values without duplicates for easy loading if required
save(betas.m.n.full, intensity.m.full, file = "/path/to/GSE105109_full_beta.RData")

# If beta/intensity values have already been generated in a previous session, they can be directly imported here.
load("/path/to/GSE105109_full_beta.RData")


## 14. Rerun ChAMP QC without duplicates, separated by assay, with CETS data
samples.m.bs <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1), "X"])
samples.m.ox <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2), "X"])

betas.m.n.full.bs <- betas.m.n.full[,samples.m.bs]
betas.m.n.full.ox <- betas.m.n.full[,samples.m.ox]

intensity.m.full.bs <- intensity.m.full[,samples.m.bs]
intensity.m.full.ox <- intensity.m.full[,samples.m.ox]

pd.m.bs <- pd.m[which(pd.m$Assay == 1), ]
pd.m.ox <- pd.m[which(pd.m$Assay == 2), ]

pd.m.bs$Slide <- as.factor(pd.m.bs$Slide) 
pd.m.ox$Slide <- as.factor(pd.m.ox$Slide) 

ssNoob.data.m.bs <- list(mset = mset.m.fm, rgSet = rgset.m, pd = pd.m.bs, intensity = intensity.m.full.bs, beta = betas.m.n.full.bs, detP = detP.m)
ssNoob.data.m.ox <- list(mset = mset.m.fm, rgSet = rgset.m, pd = pd.m.ox, intensity = intensity.m.full.ox, beta = betas.m.n.full.ox, detP = detP.m)

champ.QC(beta = ssNoob.data.m.bs$beta, pheno = pd.m.bs$post.mortem.diagnosis.ch1, resultsDir = "/path/to/GSE105109_CHAMP_QC_BS")
champ.QC(beta = ssNoob.data.m.ox$beta, pheno = pd.m.ox$post.mortem.diagnosis.ch1, resultsDir = "/path/to/GSE105109_CHAMP_QC_OX")

champ.SVD(beta = betas.m.n.full.bs, pd = pd.m.bs, resultsDir = "/path/to/GSE105109_CHAMP_SVD_BS/")
champ.SVD(beta = betas.m.n.full.ox, pd = pd.m.ox, resultsDir = "/path/to/GSE105109_CHAMP_SVD_OX/")
#Effects of glial proportion on BS data; include as covariate in models. Could also think about slide or array, although effects are less distinct.

## 15. Prep data forOxBS.MLE
#Find regions of BS and oxBS with oxbs mle

#The following code is to do a paired analysis of OxBS and BS with the ENmix package
#It requires betavalues for each sample as well as total signal (intensity)

# Select samples by assay and diagnosis
samples.m.bs.none <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Control"), "X"])
samples.m.bs.ad <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Alzheimer's disease"), "X"])

samples.m.ox.none <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Control"), "X"])
samples.m.ox.ad <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Alzheimer's disease"), "X"])

# Separate beta values and intensity values by the assay and diagnosis
betas.m.n.full.bs.none <- betas.m.n.full[,samples.m.bs.none]
betas.m.n.full.bs.ad <- betas.m.n.full[,samples.m.bs.ad]

betas.m.n.full.ox.none <- betas.m.n.full[,samples.m.ox.none]
betas.m.n.full.ox.ad <- betas.m.n.full[,samples.m.ox.ad]

intensity.m.full.bs.none <- intensity.m.full[,samples.m.bs.none]
intensity.m.full.bs.ad <- intensity.m.full[,samples.m.bs.ad]

intensity.m.full.ox.none <- intensity.m.full[,samples.m.ox.none]
intensity.m.full.ox.ad <- intensity.m.full[,samples.m.ox.ad]

# Replace column names with cases so columns match for oxBS.MLE
colnames(betas.m.n.full.bs.none) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Control"), "Case"])
colnames(betas.m.n.full.bs.ad) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Alzheimer's disease"), "Case"])

colnames(betas.m.n.full.ox.none) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Control"), "Case"])
colnames(betas.m.n.full.ox.ad) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Alzheimer's disease"), "Case"])

colnames(intensity.m.full.bs.none) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Control"), "Case"])
colnames(intensity.m.full.bs.ad) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 1 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Alzheimer's disease"), "Case"])

colnames(intensity.m.full.ox.none) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Control"), "Case"])
colnames(intensity.m.full.ox.ad) <- as.character(samples_data.m.3[which(samples_data.m.3$Assay == 2 & samples_data.m.3$post.mortem.diagnosis.ch1 == "Alzheimer's disease"), "Case"])

# Save normalized values without duplicates for easy loading if required
save(betas.m.n.full.bs.none, betas.m.n.full.bs.ad, betas.m.n.full.ox.none, betas.m.n.full.ox.ad, intensity.m.full.bs.none, intensity.m.full.bs.ad, intensity.m.full.ox.none, intensity.m.full.ox.ad, file = "/Users/JoeKochmanski/Desktop/GSE105109_RAW/GSE105109_forMLE.RData")

# If these have already been generated in a previous session, they can be directly imported here.
load("/path/to/GSE105109_forMLE.RData")


## 16. Run OxBS.MLE
# Run oxBS MLE in ENmix for each diagnosis and each tissue type
NoneMLE.m <- oxBS.MLE(beta.BS = betas.m.n.full.bs.none, beta.oxBS = betas.m.n.full.ox.none, N.BS = intensity.m.full.bs.none, N.oxBS = intensity.m.full.ox.none)
ADMLE.m <- oxBS.MLE(beta.BS = betas.m.n.full.bs.ad, beta.oxBS = betas.m.n.full.ox.ad, N.BS = intensity.m.full.bs.ad, N.oxBS = intensity.m.full.ox.ad)

# Save MLE results
save(NoneMLE.m, ADMLE.m, file = "/path/to/GSE105109_MLE_results.RData")

# Combine into matrices
dir.create("/path/to/GSE105109_RAW/oxBSMLE/")

methylcyto.m <- cbind(NoneMLE.m$`5mC`,  ADMLE.m$`5mC`)
write.table(methylcyto.m, "/path/to/GSE105109.methylcyto.m.betavals.txt", row.names = T, quote = F, sep = '\t')

hydroxymethylcyto.m <- cbind(NoneMLE.m$`5hmC`, ADMLE.m$`5hmC`)
write.table(hydroxymethylcyto.m, "/path/to/GSE105109.hydroxymethylcyto.c.betavals.txt", row.names = T, quote = F, sep = '\t')


#Save estimated 5mC and 5hmC beta values
save(methylcyto.m, hydroxymethylcyto.m, file = "/path/to/GSE105109_MLE_betas.RData")

#Load estimated 5mC and 5hmC beta values}
load("/path/to/GSE105109_MLE_betas.RData")


## 17. ChAMP QC and SVD on oxBS.MLE output

champ.QC(beta = methylcyto.m, pheno = pd.m.bs$post.mortem.diagnosis.ch1, resultsDir = "/path/to/oxBSMLE/CHAMP_QC_C_5mC")
champ.QC(beta = hydroxymethylcyto.m, pheno = pd.m.ox$post.mortem.diagnosis.ch1, resultsDir = "/path/to/oxBSMLE/CHAMP_QC_C_5hmC")

champ.SVD(beta = methylcyto.m, pd = pd.m.bs, resultsDir = "/path/to/oxBSMLE/CHAMP_SVD_C_5mC/")
champ.SVD(beta = hydroxymethylcyto.m, pd = pd.m.bs, resultsDir = "/path/to/oxBSMLE/CHAMP_SVD_C_5hmC/")
#5/29/19: Not working due to NA values (n=298) in beta value matrices.

#Replace NAs with zeroes in the data matrices:
methylcyto.m[!is.finite(methylcyto.m)] <- 0
hydroxymethylcyto.m[!is.finite(hydroxymethylcyto.m)] <- 0

champ.SVD(beta = methylcyto.m, pd = pd.m.bs, resultsDir = "/path/to/oxBSMLE/CHAMP_SVD_C_5mC/")
champ.SVD(beta = hydroxymethylcyto.m, pd = pd.m.bs, resultsDir = "/path/to/oxBSMLE/CHAMP_SVD_C_5hmC/")
#5/29/19: Not much significance for covariates in the oxBSMLE estimates

############## END OF OXBS-MLE DATA PREP AND PROCESSING #################



##############################################################################################
############### PART II: GAMLSS INTERACTION TERM MODELING ###############
##############################################################################################

##Modeling oxBS-450K data using beta regression models

#Goal: Model genome-wide changes in 5-mC and 5-hmC using mixed effects beta regression models. 

#Modeling approach:
#1. Control (n=13) vs. Braak Stage VI AD (n=22) (Only Males)

#Covariates: 
#Fixed effects = glial, age (no sex because only male) -- matches publication modeling.
#Random effect = paired ID

#After running these models, examine output to identify sites that show AD-specific differential methylation/hydroxymethylation.


##1. Filtering probes down to those with beta > 0.1 for 5-hmC and 5-mC

#Load estimated 5mC and 5hmC beta values for oxBS-450K data:
load("/path/to/oxBSMLE/GSE105109_MLE_betas.RData")

##Two beta value matrices --> methylcyto.m, hydroxymethylcyto.m
## columns are samples (n=14 control and n=22 AD) and rows are cpgs

#Install Illumina 450K array annotation file:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19", version = "3.8")

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#Illumina 450K array manifest available at: https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html

#Read in condensed version of the probe manifest that ONLY includes probe ID, chr #, start/end coordinates, and gene/regulatory annotation:
Illumina_450K <- read.csv("/path/to/HumanMethylation450_15017482_v1-2_053019.csv")

#Reformat beta values into dataframes:
hydroxymethylcyto.m.df <- as.data.frame(hydroxymethylcyto.m)
hydroxymethylcyto.m.df$probe <- rownames(hydroxymethylcyto.m.df)

methylcyto.m.df <- as.data.frame(methylcyto.m)
methylcyto.m.df$probe <- rownames(methylcyto.m.df)

#Rename sample data variables in data_frames:
colnames(hydroxymethylcyto.m.df) <- c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC",
                                      "AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC",
                                      "probe")
colnames(methylcyto.m.df) <- c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC",
                               "AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC",
                               "probe")

#Merge data frames of oxBS-MLE outputs with the 450K manifest by "probe" variable:
hydroxymethylcyto.m.df.merged <- merge(hydroxymethylcyto.m.df,Illumina_450K,by=c("probe"))
methylcyto.m.df.merged <- merge(methylcyto.m.df,Illumina_450K,by=c("probe"))
summary(methylcyto.m.df.merged)

#Create new columns of means:
hydroxymethylcyto.m.df.merged$None_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.merged[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
hydroxymethylcyto.m.df.merged$AD_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.merged[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)

methylcyto.m.df.merged$None_mean_5mC <- rowMeans(methylcyto.m.df.merged[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
methylcyto.m.df.merged$AD_mean_5mC <- rowMeans(methylcyto.m.df.merged[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)

#Filter 5-hmC and 5-mC data to remove probes with mean beta value < 0.10

library(dplyr)
hmc.m.df.merged.filter <- filter(hydroxymethylcyto.m.df.merged, None_mean_5hmC > 0.1 & AD_mean_5hmC > 0.1)
#145072 probes
mc.m.df.merged.filter <- filter(methylcyto.m.df.merged, None_mean_5mC > 0.1 & AD_mean_5mC > 0.1)
#237452 probes
summary(hmc.m.df.merged.filter)

#Now merge 5-mC and 5-hmC data for remaining:

#Join cingulate and parietal filtered probes lists by common probes
filter_bvals <- plyr::join(mc.m.df.merged.filter,hmc.m.df.merged.filter, by = "probe") #merge all 5mC and 5hmC data
row.names(filter_bvals) <- filter_bvals$probe #rename row names with probe IDs
#145072 probes in shared filtered list

summary(filter_bvals)

#Remove extraneous variables prior to modeling:
filter_bvals$probe <- NULL #remove probe variable
filter_bvals$Genome_Build <- NULL #remove Gene_symbol variable
filter_bvals$chr <- NULL 
filter_bvals$start <- NULL 
filter_bvals$UCSC_RefGene_Name <- NULL 
filter_bvals$UCSC_RefGene_Accession <- NULL 
filter_bvals$UCSC_RefGene_Group <- NULL 
filter_bvals$UCSC_CpG_Islands_Name <- NULL 
filter_bvals$Relation_to_UCSC_CpG_Island <- NULL 
filter_bvals$Phantom <- NULL 
filter_bvals$DMR <- NULL 
filter_bvals$Enhancer <- NULL 
filter_bvals$HMM_Island <- NULL 
filter_bvals$Regulatory_Feature_Name <- NULL 
filter_bvals$Regulatory_Feature_Group <- NULL 
filter_bvals$None_mean_5hmC <- NULL
filter_bvals$None_mean_5mC <- NULL 
filter_bvals$AD_mean_5hmC <- NULL 
filter_bvals$AD_mean_5mC <- NULL 

summary(filter_bvals)

filter_bvals$probe <- NULL #remove probe variable
filter_bvals$Genome_Build <- NULL #remove Gene_symbol variable
filter_bvals$chr <- NULL 
filter_bvals$start <- NULL 
filter_bvals$UCSC_RefGene_Name <- NULL 
filter_bvals$UCSC_RefGene_Accession <- NULL 
filter_bvals$UCSC_RefGene_Group <- NULL 
filter_bvals$UCSC_CpG_Islands_Name <- NULL 
filter_bvals$Relation_to_UCSC_CpG_Island <- NULL 
filter_bvals$Phantom <- NULL 
filter_bvals$DMR <- NULL 
filter_bvals$Enhancer <- NULL 
filter_bvals$HMM_Island <- NULL 
filter_bvals$Regulatory_Feature_Name <- NULL 
filter_bvals$Regulatory_Feature_Group <- NULL 
filter_bvals$None_mean_5hmC <- NULL
filter_bvals$None_mean_5mC <- NULL 
filter_bvals$AD_mean_5hmC <- NULL 
filter_bvals$AD_mean_5mC <- NULL 

summary(filter_bvals)

#Reduce data frame to complete cases:
filter_bvals.complete <- subset(filter_bvals, complete.cases(filter_bvals)) #remove probes with N/A values.
#140539 probes remaining

#create dataframe of beta values in wide format (NOT long) for downstream modeling
filter_bvals.t <- as.data.frame(t(filter_bvals.complete))

##Chunk the data frame for multiple node processing
filter_bvals.t.chnk0 <- filter_bvals.t[,c(1:5)]
filter_bvals.t.chnk1 <- filter_bvals.t[,c(1:50000)]
filter_bvals.t.chnk2 <- filter_bvals.t[,c(50001:100000)]
filter_bvals.t.chnk3 <- filter_bvals.t[,c(100001:140539)]

#Fix Beta values that are equal to zero (since zero values will mess up GAMLSS w/ beta distribution)
#THIS IS REALLY IMPORTANT FOR DOWNSTREAM TRANSFORMATIONS
if (min(filter_bvals.t.chnk0, na.rm = TRUE) <= 0)
  filter_bvals.t.chnk0[filter_bvals.t.chnk0 <= 0] <- min(filter_bvals.t.chnk0[filter_bvals.t.chnk0 > 0])
if (min(filter_bvals.t.chnk1, na.rm = TRUE) <= 0)
  filter_bvals.t.chnk1[filter_bvals.t.chnk1 <= 0] <- min(filter_bvals.t.chnk1[filter_bvals.t.chnk1 > 0])
if (min(filter_bvals.t.chnk2, na.rm = TRUE) <= 0)
  filter_bvals.t.chnk2[filter_bvals.t.chnk2 <= 0] <- min(filter_bvals.t.chnk2[filter_bvals.t.chnk2 > 0])
if (min(filter_bvals.t.chnk3, na.rm = TRUE) <= 0)
  filter_bvals.t.chnk3[filter_bvals.t.chnk3 <= 0] <- min(filter_bvals.t.chnk3[filter_bvals.t.chnk3 > 0])

#################################################################
##2. Establish Meta Data Modeling Variables

#Read in meta data that also includes the glial proportions
metadata <- read.csv("/path/to/GSE105109.pheno.male.filtered.v4.053019.csv", header = T, check.names = F, stringsAsFactors = F)

#Get the glial cell composition
glial <- metadata$glial
age <- metadata$age.at.death.ch1

#Build the design
Disease_cat <- as.factor(c(rep(1, 13), rep(2, 22),rep(1, 13), rep(2, 22)))
Disease_cat 
length(Disease_cat)
DNA_mod_cat <- as.factor(c(rep(1, 35), rep(2, 35)))
length(DNA_mod_cat)

#Paired ID as a random effect (Assay = A or B, depending on 5-mC or 5-hmC)
randeff <- as.factor(c(1:35,1:35))
randeff
length(randeff)

###########################################################
######### Beta regression modeling functions ##############
###########################################################

#Load required packages
library(plyr)
library(dplyr)
library(gamlss)
library(parallel)

#NOTE: For this filtered data, data is no longer zero-inflated. As a result, use different family for modeling function!
#Example: beta binomial distribution


#############################################################
######### Model: Control vs. AD Comparison ##############

#Beta regression with a logit link function
Ilm450k.beta.fit.full <- function(x, glial, age, randeff) {
  x <- as.data.frame(x)
  x$randeff <- randeff
  x$glial <- glial
  x$age <- age
  
  #Define data
  data <- x[1:70,]
  
  #Fit models
  #Fit with beta regression with C_None as the reference
  none.fit <- gamlss(x ~ Disease_cat*DNA_mod_cat + glial + age, random = ~1 | randeff, data = data, family = BE, trace = F)
  
  #Get the beta estimates for fold-change calcs
  none.fit.sum <- as.data.frame(summary(none.fit))[6,c(1:4)] #Save estimate, se, t-value, and P-value [1:4] for beta coefficient 6 (interaction term)
  #Build the output df
  #cnone.fit.sum$pval <- rbind(cnone.lrt.limbic, cnone.lrt.neocortical)
  
  return(none.fit.sum)
}

#Test the model on small chunk of probes
methylcyto.filter.fit.0 <- mclapply(as.data.frame(filter_bvals.t.chnk0), function(x) Ilm450k.beta.fit.full(x, glial, age, randeff), mc.preschedule = F, mc.cores = 4) #TEST! (5/30/19)
summary(methylcyto.filter.fit.0)
methylcyto.filter.fit.0$cg00000109 #Example for how to view model coefficients for an individual probe in a specific tissue type
#Try making data frame out of individual CpG probe in list:
methylcyto.filter.fit.0.df3 <- ldply (methylcyto.filter.fit.0, data.frame) #Note: this removes the rownames; need to recreate as separate variable

# Scaling up to a large chunk of CpG sites (n=50,000):
methylcyto.filter.fit.1 <- mclapply(as.data.frame(filter_bvals.t.chnk1), function(x) Ilm450k.beta.fit.full(x, glial, age, randeff), mc.preschedule = F, mc.cores = 4) 
methylcyto.filter.fit.2 <- mclapply(as.data.frame(filter_bvals.t.chnk2), function(x) Ilm450k.beta.fit.full(x, glial, age, randeff), mc.preschedule = F, mc.cores = 4)
methylcyto.filter.fit.3 <- mclapply(as.data.frame(filter_bvals.t.chnk3), function(x) Ilm450k.beta.fit.full(x, glial, age, randeff), mc.preschedule = F, mc.cores = 4) 

summary(methylcyto.filter.fit.1)
summary(methylcyto.filter.fit.2)
summary(methylcyto.filter.fit.3)

#Now create long format data frames of results for all probes:
Variable <- c("Disease_cat2.DNA_mod_cat2")
Variable_names <- c(rep(Variable, 50000)) #Create replicated variable name vector
Variable_names2 <- c(rep(Variable, 40539)) #Create replicated variable name vector

methylcyto.filter.fit.1.df3 <- ldply (methylcyto.filter.fit.1, data.frame) #Note: this removes the rownames; need to recreate as separate variable
methylcyto.filter.fit.1.df3$variables <- Variable_names #Add variable name vector to data frame of results

methylcyto.filter.fit.2.df3 <- ldply (methylcyto.filter.fit.2, data.frame) #Note: this removes the rownames; need to recreate as separate variable
methylcyto.filter.fit.2.df3$variables <- Variable_names #Add variable name vector to data frame of results

methylcyto.filter.fit.3.df3 <- ldply (methylcyto.filter.fit.3, data.frame) #Note: this removes the rownames; need to recreate as separate variable
methylcyto.filter.fit.3.df3$variables <- Variable_names2 #Add variable name vector to data frame of results

##Now bind the separate data frames together with rbind:
methylcyto.filter.fit.full.df1 <- rbind(methylcyto.filter.fit.1.df3,methylcyto.filter.fit.2.df3)
methylcyto.filter.fit.full.df <- rbind(methylcyto.filter.fit.full.df1,methylcyto.filter.fit.3.df3) #Final dataframe of results
summary(methylcyto.filter.fit.full.df)

#Given a set of p-values, returns p-values adjusted using FDR method in p.adjust function
p1 <- methylcyto.filter.fit.full.df$Pr...t..

fdr <- p.adjust(p1, method = "BH", n = length(p1)) #Benjamini-Hochberg FDR adjustment for parietal p-values

#Append the FDR values to the dataset.
methylcyto.filter.fit.full.df$fdr <- fdr

#Filter down to probes with standard error less than 10 to remove the subpopulation of probes with standard errors in the 10000s. These are poorly modeled in this analysis.
methylcyto.filter.fit.full.df1 <- methylcyto.filter.fit.full.df[(methylcyto.filter.fit.full.df$Std..Error <10), ] 
# 140539 probes

#Filter by p-value (p<0.001):
methylcyto.filter.fit.full.df.pval <- methylcyto.filter.fit.full.df1[(methylcyto.filter.fit.full.df1$Pr...t.. < 0.001), ] 
# 5955 probes

#Filter by FDR value (p<0.05):
methylcyto.filter.fit.full.df.fdr <- methylcyto.filter.fit.full.df1[(methylcyto.filter.fit.full.df1$fdr < 0.05), ] 
methylcyto.filter.fit.full.df.fdr$probe <- methylcyto.filter.fit.full.df.fdr$.id
# 14183 probes

##Now examine data across the significant sites:

#Before I can visualize the data at these twelve sites, I need to annotate the probe IDs.

#Install Illumina 450K array annotation file:
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#Illumina 450K array manifest available at: https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html

#Read in condensed version of the probe manifest that ONLY includes probe ID, chr #, start/end coordinates, and gene/regulatory annotation:
Illumina_450K <- read.csv("/path/to/HumanMethylation450_15017482_v1-2_053019.csv")

#Establish probe ID as rowname for each dataset:
methylcyto.filter.fit.full.df.fdr$probe <- methylcyto.filter.fit.full.df.fdr$.id

#Merge data frames of gamlss interaction term model outputs with the EPIC manifest by "probe" variable:
methylcyto.filter.fit.full.df.fdr.merged <- merge(methylcyto.filter.fit.full.df.fdr,Illumina_450K,by=c("probe"))

#Now remove the extra columns introduced from EPIC manifest
#methylcyto.filter.fit.full.df.fdr.merged <- methylcyto.filter.fit.full.df.fdr.merged[ -c(19:119) ]

#Sort the dataframe by FDR values:
methylcyto.filter.fit.full.df.fdr.merged.sorted <- methylcyto.filter.fit.full.df.fdr.merged[order(methylcyto.filter.fit.full.df.fdr.merged$fdr),] 

#Quickly check raw beta values to see whether significance is reflected in values.

#Pull out raw beta values for significant probes (by FDR value)
sig.sites1 <- methylcyto.filter.fit.full.df.fdr.merged.sorted$probe

#subset the larger dataset by sig.sites variable to view data for ALL significant probes:
hydroxymethylcyto.m.df.sig.sites1 <- subset(hydroxymethylcyto.m.df, (probe %in% sig.sites1))
methylcyto.m.df.sig.sites1 <- subset(methylcyto.m.df, (probe %in% sig.sites1))

#Create new columns of means for filtered probes:
hydroxymethylcyto.m.df.sig.sites1$None_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites1[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
hydroxymethylcyto.m.df.sig.sites1$AD_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites1[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)

methylcyto.m.df.sig.sites1$None_mean_5mC <- rowMeans(methylcyto.m.df.sig.sites1[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
methylcyto.m.df.sig.sites1$AD_mean_5mC <- rowMeans(methylcyto.m.df.sig.sites1[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)

#Add these new columns of means to the spreadsheet of filtered/significant probes:
methylcyto.filter.fit.full.df.fdr.merged.sorted$None_mean_beta_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites1[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
methylcyto.filter.fit.full.df.fdr.merged.sorted$AD_mean_beta_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites1[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)
methylcyto.filter.fit.full.df.fdr.merged.sorted$None_mean_beta_5mC <- rowMeans(methylcyto.m.df.sig.sites1[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
methylcyto.filter.fit.full.df.fdr.merged.sorted$AD_mean_beta_5mC <- rowMeans(methylcyto.m.df.sig.sites1[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)

#Write .csv files of the results:
write.csv(methylcyto.filter.fit.full.df.fdr.merged.sorted, "/path/to/GSE105109_oxBS-450K_controlvsAD_5mC.5hmC.interaction.results.nopmi.fullFDR_053019.csv")

##############################################################################################
############ PART III: GAMLSS SEPARATE 5-mC and 5-hmC BETA REGRESSION MODELING ###############
##############################################################################################


##Next steps:
# 1. Run separate beta regression models (w/o interaction term) for 5-mC and 5-hmC.

#Read in meta data that also includes the glial proportions
metadata <- read.csv("/path/to/GSE105109.pheno.male.filtered.v4.053019.csv", header = T, check.names = F, stringsAsFactors = F)

#Get the glial cell composition and age variables
glial <- metadata$glial
age <- metadata$age.at.death.ch1

#Build the design
Disease_cat.sep <- as.factor(c(rep(1, 13), rep(2, 22)))
Disease_cat.sep 
length(Disease_cat.sep)

#Load required packages
library(plyr)
library(dplyr)
library(gamlss)
library(parallel)

####################################################################
######### Model: Control vs. AD Comparison - 5mC ONLY ##############

#Beta regression with a logit link function
Ilm450k.beta.fit.sep <- function(x, glial, age) {
  x <- as.data.frame(x)
  x$glial <- glial
  x$age <- age
  
  #Define data
  data.5mC <- x[1:35,]
  data.5hmC <- x[36:70,]
  
  #Fit models
  #Fit with beta regression with C_None as the reference
  mc.fit <- gamlss(x ~ Disease_cat.sep + glial + age, data = data.5mC, family = BE, trace = F)
  
  #Get the beta estimates for fold-change calcs
  mc.fit.sum <- as.data.frame(summary(mc.fit))[2,c(1:4)] #Save estimate, se, t-value, and P-value [1:4] for beta coefficient 1 (disease term)
  #Build the output df
  #cnone.fit.sum$pval <- rbind(cnone.lrt.limbic, cnone.lrt.neocortical)
  
  #Fit with beta regression with C_None as the reference
  hmc.fit <- gamlss(x ~ Disease_cat.sep + glial + age, data = data.5hmC, family = BE, trace = F)
  
  #Get the beta estimates for fold-change calcs
  hmc.fit.sum <- as.data.frame(summary(hmc.fit))[2,c(1:4)] #Save estimate, se, t-value, and P-value [1:4] for beta coefficient 1 (disease term)
  
  return(list(model.5mC = mc.fit.sum,
              model.5hmC = hmc.fit.sum))
}

#Test the model on small chunk of probes
methylcyto.filter.fit.0 <- mclapply(as.data.frame(filter_bvals.t.chnk0), function(x) Ilm450k.beta.fit.sep(x, glial, age), mc.preschedule = F, mc.cores = 4) #TEST! (5/30/19)
summary(methylcyto.filter.fit.0)
methylcyto.filter.fit.0$cg00000109$model.5hmC #Example for how to view model coefficients for an individual probe in a specific tissue type
methylcyto.filter.fit.0$cg00000109$model.5mC #Example for how to view model coefficients for an individual probe in a specific tissue type
#Try making data frame out of individual CpG probe in list:
methylcyto.filter.fit.0.df3 <- ldply (methylcyto.filter.fit.0, data.frame) #Note: this removes the rownames; need to recreate as separate variable
#Establish probe ID as rowname for each dataset:
#methylcyto.filter.fit.0.df3$probe <- methylcyto.filter.fit.0.df3$.id
#Merge data frames of gamlss separate term model outputs with the EPIC manifest by "probe" variable:
#methylcyto.filter.fit.0.df3.merged <- merge(methylcyto.filter.fit.0.df3,Illumina_450K,by=c("probe"))
#Establish probe IDs to pull out from raw beta data:
#sig.sites.0 <- methylcyto.filter.fit.0.df3.merged$probe
#subset the larger beta dataset by sig.sites variable to view data for ALL significant probes:
#methylcyto.m.df.sig.sites0 <- subset(methylcyto.m.df, (probe %in% sig.sites.0))
#hydroxymethylcyto.m.df.sig.sites0 <- subset(hydroxymethylcyto.m.df, (probe %in% sig.sites.0))
#Create new columns of means for filtered probes:
#methylcyto.filter.fit.0.df3.merged$None_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites0[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
#methylcyto.filter.fit.0.df3.merged$AD_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites0[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)
#methylcyto.filter.fit.0.df3.merged$None_mean_5mC <- rowMeans(methylcyto.m.df.sig.sites0[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
#methylcyto.filter.fit.0.df3.merged$AD_mean_5mC <- rowMeans(methylcyto.m.df.sig.sites0[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)

#Add these new columns of means to the spreadsheet of filtered/significant probes:
hmc.filter.fit.sep.full.df.pval.merged.sorted$None_mean_beta_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites2[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
hmc.filter.fit.sep.full.df.pval.merged.sorted$AD_mean_beta_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites2[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)
mc.filter.fit.sep.full.df.pval.merged.sorted$None_mean_beta_5mC <- rowMeans(methylcyto.m.df.sig.sites2[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
mc.filter.fit.sep.full.df.pval.merged.sorted$AD_mean_beta_5mC <- rowMeans(methylcyto.m.df.sig.sites2[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)


#5/30/19 - Now try scaling up to a large chunk of CpG sites (n=50,000):
methylcyto.filter.fit.sep.1 <- mclapply(as.data.frame(filter_bvals.t.chnk1), function(x) Ilm450k.beta.fit.sep(x, glial, age), mc.preschedule = F, mc.cores = 4) 
methylcyto.filter.fit.sep.2 <- mclapply(as.data.frame(filter_bvals.t.chnk2), function(x) Ilm450k.beta.fit.sep(x, glial, age), mc.preschedule = F, mc.cores = 4)
methylcyto.filter.fit.sep.3 <- mclapply(as.data.frame(filter_bvals.t.chnk3), function(x) Ilm450k.beta.fit.sep(x, glial, age), mc.preschedule = F, mc.cores = 4) 

summary(methylcyto.filter.fit.sep.1)
summary(methylcyto.filter.fit.sep.2)
summary(methylcyto.filter.fit.sep.3)

#Now create long format data frames of results for all probes:
Variable <- c("Disease_cat.sep2")
Variable_names <- c(rep(Variable, 50000)) #Create replicated variable name vector
Variable_names2 <- c(rep(Variable, 40539)) #Create replicated variable name vector

methylcyto.filter.fit.sep.1.df3 <- ldply (methylcyto.filter.fit.sep.1, data.frame) #Note: this removes the rownames; need to recreate as separate variable
methylcyto.filter.fit.sep.1.df3$variables <- Variable_names #Add variable name vector to data frame of results

methylcyto.filter.fit.sep.2.df3 <- ldply (methylcyto.filter.fit.sep.2, data.frame) #Note: this removes the rownames; need to recreate as separate variable
methylcyto.filter.fit.sep.2.df3$variables <- Variable_names #Add variable name vector to data frame of results

methylcyto.filter.fit.sep.3.df3 <- ldply (methylcyto.filter.fit.sep.3, data.frame) #Note: this removes the rownames; need to recreate as separate variable
methylcyto.filter.fit.sep.3.df3$variables <- Variable_names2 #Add variable name vector to data frame of results

##Now bind the separate data frames together with rbind:
methylcyto.filter.fit.sep.full.df1 <- rbind(methylcyto.filter.fit.sep.1.df3,methylcyto.filter.fit.sep.2.df3)
methylcyto.filter.fit.sep.full.df <- rbind(methylcyto.filter.fit.sep.full.df1,methylcyto.filter.fit.sep.3.df3) #Final dataframe of results
summary(methylcyto.filter.fit.sep.full.df)


#Pull out p-values for 5mC and 5hmC modeling:
p2 <- methylcyto.filter.fit.sep.full.df$model.5mC.Pr...t..
p3 <- methylcyto.filter.fit.sep.full.df$model.5hmC.Pr...t..
p4 <- c(p2,p3) #Combine the p-values into a shared vector (complete list)
length(p4)

#Given a set of p-values, p.adjust() functions returns p-values adjusted using correction method (e.g. Benjamini-Hochberg FDR)
fdr.5mC <- p.adjust(p2, method = "BH", n = length(p2)) #Benjamini-Hochberg FDR adjustment for 5mC p-values
fdr.5hmC <- p.adjust(p3, method = "BH", n = length(p3)) #Benjamini-Hochberg FDR adjustment for 5hmC p-values

#Append the FDR values to the dataset.
methylcyto.filter.fit.sep.full.df$fdr.5mC <- fdr.5mC
methylcyto.filter.fit.sep.full.df$fdr.5hmC <- fdr.5hmC

#Filter down to probes with standard error less than 10 to remove the subpopulation of probes with standard errors in the 10000s. These are poorly modeled in this analysis.
methylcyto.filter.fit.sep.full.df1 <- methylcyto.filter.fit.sep.full.df[(methylcyto.filter.fit.sep.full.df$model.5mC.Std..Error <10 & methylcyto.filter.fit.sep.full.df$model.5hmC.Std..Error <10), ] 
# 140539 probes

#Filter by p-value (p<0.001):
mc.filter.fit.sep.full.df.pval <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$model.5mC.Pr...t.. < 0.001), ] 
# 232 probes
hmc.filter.fit.sep.full.df.pval <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$model.5hmC.Pr...t.. < 0.001), ] 
# 568 probes

#Filter by FDR value (p<0.05):
mc.filter.fit.sep.full.df.fdr <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$fdr.5mC < 0.05), ] 
#2 probes
hmc.filter.fit.sep.full.df.fdr <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$fdr.5hmC < 0.05), ] 
#0 probes

#Filter by FDR value (p<0.10):
mc.filter.fit.sep.full.df.fdr <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$fdr.5mC < 0.10), ] 
#2 probes
hmc.filter.fit.sep.full.df.fdr <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$fdr.5hmC < 0.10), ] 
#1 probe


#Re-establish "probe" variable from .id variable
mc.filter.fit.sep.full.df.fdr$probe <- mc.filter.fit.sep.full.df.fdr$.id
hmc.filter.fit.sep.full.df.fdr$probe <- hmc.filter.fit.sep.full.df.fdr$.id
# 2 probes

##Now examine data across the significant sites:

#Before I can visualize the data at these twelve sites, I need to annotate the probe IDs.

#Install Illumina 450K array annotation file:
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

#Illumina 450K array manifest available at: https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html

#Read in condensed version of the probe manifest that ONLY includes probe ID, chr #, start/end coordinates, and gene/regulatory annotation:
Illumina_450K <- read.csv("/path/to/HumanMethylation450_15017482_v1-2_053019.csv")

#Establish probe ID as rowname for each dataset:
mc.filter.fit.sep.full.df.pval$probe <- mc.filter.fit.sep.full.df.pval$.id
hmc.filter.fit.sep.full.df.pval$probe <- hmc.filter.fit.sep.full.df.pval$.id

#Merge data frames of gamlss separate term model outputs with the EPIC manifest by "probe" variable:
mc.filter.fit.sep.full.df.pval.merged <- merge(mc.filter.fit.sep.full.df.pval,Illumina_450K,by=c("probe"))
hmc.filter.fit.sep.full.df.pval.merged <- merge(hmc.filter.fit.sep.full.df.pval,Illumina_450K,by=c("probe"))

#Sort the dataframe by FDR values:
mc.filter.fit.sep.full.df.pval.merged.sorted <- mc.filter.fit.sep.full.df.pval.merged[order(mc.filter.fit.sep.full.df.pval.merged$model.5mC.Pr...t..),] 
hmc.filter.fit.sep.full.df.pval.merged.sorted <- hmc.filter.fit.sep.full.df.pval.merged[order(hmc.filter.fit.sep.full.df.pval.merged$model.5hmC.Pr...t..),] 

#Could also perform these steps for FDR (but don't bother since it's only a couple probes...)
#Merge data frames of gamlss separate term model outputs with the EPIC manifest by "probe" variable:
#mc.filter.fit.sep.full.df.fdr.merged <- merge(mc.filter.fit.sep.full.df.fdr,Illumina_450K,by=c("probe"))
#hmc.filter.fit.sep.full.df.fdr.merged <- merge(hmc.filter.fit.sep.full.df.fdr,Illumina_450K,by=c("probe"))

#Sort the dataframe by FDR values:
#mc.filter.fit.full.df.fdr.merged.sorted <- mc.filter.fit.sep.full.df.fdr.merged[order(mc.filter.fit.sep.full.df.fdr.merged$fdr.5mC),] 
#hmc.filter.fit.full.df.fdr.merged.sorted <- hmc.filter.fit.sep.full.df.fdr.merged[order(hmc.filter.fit.sep.full.df.fdr.merged$fdr.5mC),] 


#Quickly check raw beta values to see whether significance is reflected in values.
#Pull out raw beta values for significant probes (by p-value < 0.001 due to lack of results by FDR cutoff)

#Establish probe IDs to pull out from raw beta data:
sig.sites.mc <- mc.filter.fit.sep.full.df.pval.merged.sorted$probe
sig.sites.hmc <- hmc.filter.fit.sep.full.df.pval.merged.sorted$probe

#subset the larger beta dataset by sig.sites variable to view data for ALL significant probes:
hydroxymethylcyto.m.df.sig.sites2 <- subset(hydroxymethylcyto.m.df, (probe %in% sig.sites.hmc))
methylcyto.m.df.sig.sites2 <- subset(methylcyto.m.df, (probe %in% sig.sites.mc))

#Create new columns of means for filtered probes:
hydroxymethylcyto.m.df.sig.sites2$None_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites2[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
hydroxymethylcyto.m.df.sig.sites2$AD_mean_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites2[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)

methylcyto.m.df.sig.sites2$None_mean_5mC <- rowMeans(methylcyto.m.df.sig.sites2[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
methylcyto.m.df.sig.sites2$AD_mean_5mC <- rowMeans(methylcyto.m.df.sig.sites2[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)

#Add these new columns of means to the spreadsheet of filtered/significant probes:
hmc.filter.fit.sep.full.df.pval.merged.sorted$None_mean_beta_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites2[c("None1_5hmC","None2_5hmC","None3_5hmC","None4_5hmC","None5_5hmC","None6_5hmC","None7_5hmC","None8_5hmC","None9_5hmC","None10_5hmC","None11_5hmC","None12_5hmC","None13_5hmC")], na.rm=TRUE)
hmc.filter.fit.sep.full.df.pval.merged.sorted$AD_mean_beta_5hmC <- rowMeans(hydroxymethylcyto.m.df.sig.sites2[c("AD1_5hmC","AD2_5hmC","AD3_5hmC", "AD4_5hmC","AD5_5hmC","AD6_5hmC","AD7_5hmC","AD8_5hmC","AD9_5hmC","AD10_5hmC","AD11_5hmC","AD12_5hmC","AD13_5hmC","AD14_5hmC","AD15_5hmC","AD16_5hmC","AD17_5hmC","AD18_5hmC","AD19_5hmC","AD20_5hmC","AD21_5hmC","AD22_5hmC")], na.rm=TRUE)
mc.filter.fit.sep.full.df.pval.merged.sorted$None_mean_beta_5mC <- rowMeans(methylcyto.m.df.sig.sites2[c("None1_5mC","None2_5mC","None3_5mC","None4_5mC","None5_5mC","None6_5mC","None7_5mC","None8_5mC","None9_5mC","None10_5mC","None11_5mC","None12_5mC","None13_5mC")], na.rm=TRUE)
mc.filter.fit.sep.full.df.pval.merged.sorted$AD_mean_beta_5mC <- rowMeans(methylcyto.m.df.sig.sites2[c("AD1_5mC","AD2_5mC","AD3_5mC", "AD4_5mC","AD5_5mC","AD6_5mC","AD7_5mC","AD8_5mC","AD9_5mC","AD10_5mC","AD11_5mC","AD12_5mC","AD13_5mC","AD14_5mC","AD15_5mC","AD16_5mC","AD17_5mC","AD18_5mC","AD19_5mC","AD20_5mC","AD21_5mC","AD22_5mC")], na.rm=TRUE)

#Write .csv files of the results:
write.csv(hmc.filter.fit.sep.full.df.pval.merged.sorted, "/Users/JoeKochmanski/Desktop/GSE105109_RAW/GSE105109_oxBS-450K_controlvsAD_ONLY5hmC.results.nopmi.pval_053119.csv")
write.csv(mc.filter.fit.sep.full.df.pval.merged.sorted, "/Users/JoeKochmanski/Desktop/GSE105109_RAW/GSE105109_oxBS-450K_controlvsAD_ONLY5mC.results.nopmi.pval_053119.csv")


# 2. Examine the overlap of those two separate analyses.

#Pull out sites where BOTH 5mC and 5hmC have FDR < 0.10:
mc.filter.fit.sep.full.df.fdr.full <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$fdr.5mC < 0.10 & methylcyto.filter.fit.sep.full.df1$fdr.5hmC < 0.10), ] 
#0 probes

#Pull out sites where BOTH 5mC and 5hmC have p-value < 0.001:
mc.filter.fit.sep.full.df.pval.full <- methylcyto.filter.fit.sep.full.df1[(methylcyto.filter.fit.sep.full.df1$model.5mC.Pr...t.. < 0.001 & methylcyto.filter.fit.sep.full.df1$model.5hmC.Pr...t.. < 0.001), ] 
#68 probes

#Write csv of the overlapping probes:
write.csv(mc.filter.fit.sep.full.df.pval.full, "/path/to/GSE105109_oxBS-450K_controlvsAD_BOTH5mC+5hmC.results.nopmi.pval_053119.csv")

#Re-establish "probe" variable from .id variable
mc.filter.fit.sep.full.df.pval.full$probe <- mc.filter.fit.sep.full.df.pval.full$.id

#Establish probe IDs to pull out from raw beta data:
sig.sites.full <- mc.filter.fit.sep.full.df.pval.full$probe
# probes


# 3. Check how many of the overlapping sites show up (or don't) in our paired interaction term models

#Check for overlap between the interaction term models and the probe IDs that overlapped in separate models:

sig.sites.overlap <- merge(mc.filter.fit.sep.full.df.pval.full, methylcyto.filter.fit.full.df.fdr.merged.sorted, by = "probe")
#All 68 sites overlap! Missing out on a HUGE number of potentially interesting sites by analyzing separately! 

write.csv(sig.sites.overlap, "/path/to/GSE105109_oxBS-450K_controlvsAD_5mC+5hmCoverlap_SepModels.results.nopmi.pval_053119.csv")


### GOAL: Check overlap between separate modeling approach probes (that DON'T overlap between 5-mC and 5-hmC) and interaction term probes!

#Calculate the number of unique probe IDs with significant probes in parallel, separated modeling (NOT overlapping between separate models)
#as well as the paired interaction term modeling:

# 4. Read in datasets for separate modeling approach and interaction term modeling approach:

#Significant differentially methylated probes:
mc.m.df <- read.csv("/path/to/GSE105109_oxBS-450K_controlvsAD_ONLY5mC.results.nopmi.pval_053119.csv")

#Significant differentially hydroxymethylated probes:
hmc.m.df <- read.csv("/path/to/GSE105109_oxBS-450K_controlvsAD_ONLY5hmC.results.nopmi.pval_053119.csv")

#Significant interaction term hhydroxymethylated AND methylated probes:
hmc.mc.m.int.df <- read.csv("/path/to/GSE105109_oxBS-450K_controlvsAD_5mC.5hmC.interaction.results.nopmi.fullFDR_053019.csv")

# 5. Examine the overlap of the two separate analyses to establish probe IDs that we DON'T want to compare.

#Pull out sites where BOTH 5mC and 5hmC have p-value < 0.001:
mc.m.df.pval.full <- mc.m.df[(mc.m.df$model.5mC.Pr...t.. < 0.001 & mc.m.df$model.5hmC.Pr...t.. < 0.001), ] 
#68 probes

568-68
232-68

probe.remove <- mc.m.df.pval.full$probe
class(probe.remove)

# 6. Check how many of the overlapping sites show up (or don't) in our paired interaction term models

#Check for overlap between the interaction term models and the probe IDs that overlapped in separate models:

sig.sites.overlap.mc <- merge(mc.m.df, hmc.mc.m.int.df, by = "probe") #195 probes
sig.sites.overlap.hmc <- merge(hmc.m.df, hmc.mc.m.int.df, by = "probe") #551 probes

#Now remove rows from these data frames where the probe ID = remove probe list (significant overlapping separate model probes AND interaction probes)

sig.sites.overlap.mc2 <- sig.sites.overlap.mc[ ! sig.sites.overlap.mc$probe %in% probe.remove, ] #127 probes (out of 164 probes)
sig.sites.overlap.hmc2 <- sig.sites.overlap.hmc[ ! sig.sites.overlap.hmc$probe %in% probe.remove, ] #483 probes (out of 500 probes)
127+483 #610 probes overlap between interaction term modeling and ONLY ONE of the separate models.
(127+483)/(164+500) #[1] 0.9186747 --> 91.9% of probes that were significant in ONLY ONE model were also significant in interaction term modeling.
14183-127-483-68 #13505 probes would not have been identified with separate models alone!

#Number of probes missed by interaction term modeling
#5-mC:
164-127 #37 probes
#5-hmC:
500-483 #17 probes

#Total:
37+17 #54 probes missed by interaction term modeling.
