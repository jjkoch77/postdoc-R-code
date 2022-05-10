#5/5/21 - Code for viewing existing sleuth results saved as .RData files

#First, you will need to load the R packages.
install.packages("remotes")
install.packages("WGCNA")
source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager")
BiocManager::install("remotes")
BiocManager::install(c("pachterlab/sleuth"))
BiocManager::install("WGCNA")
BiocManager::install("tidyverse")
library(WGCNA)
library(sleuth)
library(tidyverse)

# Next, load the previously saved RData files of RNA-seq sleuth results. Make sure to change /path/to/directory/ to correspond to the location of the
# files on your local computer.

#Log-likelihood ratio test results (differential expression hypothesis testing):
load("/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/AnnaStoll_DGE_sleuth_object_female_LRT_040921.RData")
load("/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/AnnaStoll_DGE_sleuth_object_male_LRT_040921.RData")
load("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/AnnaStoll_DGE_sleuth_object_male_LRT_040921.RData")
#Wald test results (only used for extracting beta coefficients)
load("/path/to/directory/P3_RNA-seq/AnnaStoll_DGE_sleuth_object_female_WT_040921.RData")
load("/path/to/directory/P3_RNA-seq/AnnaStoll_DGE_sleuth_object_male_WT_040921.RData")

#Now that data is loaded, the LRT expression testing results can be viewed using sleuth_live, a shiny app built into the sleuth package.
sleuth_live(so.F) #Female data
sleuth_live(so.M) #Male data

#Note: Male and female data were saved separately because all tests were stratified by sex. Can only open one dataset
# at a time in sleuth_live shiny app.

#Turn off sleuth_live interface when done.
dev.off()

###################### WGCNA Analysis #######################

#Convert sleuth to tpm table for downstream processing
#Example code:
#sleuth_to_matrix(obj, which_df, which_units)
#parameters:
#obj - a sleuth object
#which_df - character vector of length one. Which type of data to use ("obs_norm" or "obs_raw")
#which_units - character vector of length one. Which units to use ("tpm" or "est_counts")

so.F.tpm <- sleuth_to_matrix(so.F, which_df="obs_norm", which_units="tpm")
so.M.tpm <- sleuth_to_matrix(so.M, which_df="obs_norm", which_units="tpm")

colnames(so.M.tpm)
head(so.M.tpm)

#Note: We suggest removing features whose counts are consistently low (for example, removing all features that have a count of less than say 10 in more than 90% of the samples) 
#because such low-expressed features tend to reflect noise and correlations based on counts that are mostly zero aren't really meaningful. 
#The actual thresholds should be based on experimental design, sequencing depth and sample counts.

#Establish data frame of results:
so.M.tpm1 <- as.data.frame(so.M.tpm)

#As one option, we can filter out all transcripts where less than 50% have < 1 tpm:
so.M.tpm2 <- so.M.tpm1[rowSums(so.M.tpm1 >= 1) >= (ncol(so.M.tpm1)/2), ] 
#35527 transcripts

#Alternatively, we can be more stringent, and only allow samples where all (except M9C, which is a known outlier) have counts > 1 tpm:
so.M.tpm2 <- so.M.tpm1[which(so.M.tpm1$sampM10C_quant > 1 & so.M.tpm1$sampM11C_quant > 1 & so.M.tpm1$sampM1P_quant > 1 & so.M.tpm1$sampM2P_quant > 1 & 
                               so.M.tpm1$sampM3P_quant > 1 & so.M.tpm1$sampM5P_quant > 1 & so.M.tpm1$sampM8C_quant > 1), ]
#27564 (out of 113100 total obs) transcripts have > 1 tpm for all samples.


#Now follow the WGCNA tutorial
#Tutorial available here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Display the current working directory
#getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir); 
# Load the package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
maleData = so.M.tpm2;
maleData = as.data.frame(maleData)
# Take a quick look at what is in the data sets (caution, longish output):
dim(maleData)
names(maleData)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datExpr0 = as.data.frame(t(maleData));
names(datExpr0) = rownames(maleData);
rownames(datExpr0) = names(maleData);

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 100000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 100000, minSize = 4)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================


traitData = read.csv("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/RNAseq_Sortwell_Meta_WGCNA_051021.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData
allTraits = allTraits[, c(1,2,5)];
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
maleSamples = rownames(datExpr);
traitRows = match(maleSamples, allTraits$sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];

collectGarbage();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================


save(datExpr, datTraits, file = "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/MaleBrain-01-dataInput.RData")

#=====================================================================================
#
#  Code chunk 10 (If necessary)
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/MaleBrain-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames


#=====================================================================================
#
#  Code chunk 11
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#power = 9 appears to provide > 0.9 scale free topology fit for the male data.

#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "maleMouseTOM", 
                       verbose = 3)


#=====================================================================================
#
#  Code chunk 13
#
#=====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk 14
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/MaleBrain-02-networkConstruction-auto.RData")


#=====================================================================================
#
#  Code chunk 15 (if necessary)
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); 
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/MaleBrain-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/MaleBrain-02-networkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  Code chunk 16
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#Order modules by relationship (p-value) with treatment variable
moduleTraitPvalue1 <- as.data.frame(moduleTraitPvalue)
moduleTraitPvalue1 <- moduleTraitPvalue1[order(moduleTraitPvalue1$treatment),]

#Now filter down to those modules with p-value < 0.05 for treatment variable:
moduleTraitPvalue2 <- moduleTraitPvalue1[which(moduleTraitPvalue1$treatment<0.05),]
#n=5 out of 85 modules with p-value <0.05:
#MEturquoise - p-value = 0.0006717653
#MEbrown - p-value = 0.0145116281
#MEthistle - p-value = 0.0254840687
#MEpurple - p-value = 0.0409220157
#MEskyblue - p-value = 0.0442125436

#The turquoise module is the most significant module! Then brown, thistle, purple, and skyblue, respectively.

#=====================================================================================
#
#  Code chunk 17
#
#=====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#=====================================================================================
#
#  Code chunk 18
#
#=====================================================================================


# Define variable treatment containing the treatment column of datTrait
treatment = as.data.frame(datTraits$treatment);
names(treatment) = "treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, treatment, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(treatment), sep="");
names(GSPvalue) = paste("p.GS.", names(treatment), sep="");


#=====================================================================================
#
#  Code chunk 19
#
#=====================================================================================
#Male data: turquoise module is the most significant module! Then brown, thistle, purple, and skyblue, respectively.


module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PFF treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PFF treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "thistle"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PFF treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "purple"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PFF treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

module = "skyblue"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for PFF treatment",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#=====================================================================================
#
#  Code chunk 20
#
#=====================================================================================


names(datExpr)


#=====================================================================================
#
#  Code chunk 21
#
#=====================================================================================
#Male data: turquoise module is the most significant module! Then brown, thistle, purple, and skyblue, respectively.


names(datExpr)[moduleColors=="turquoise"]
names(datExpr)[moduleColors=="brown"]
names(datExpr)[moduleColors=="thistle"]
names(datExpr)[moduleColors=="purple"]
names(datExpr)[moduleColors=="skyblue"]


#=====================================================================================
#
#  Code chunk 22
#
#=====================================================================================


#annot = read.csv(file = "GeneAnnotation.csv");
BiocManager::install("rtracklayer")
library(rtracklayer)
gtf <- rtracklayer::import('/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/Rattus_norvegicus.Rnor_6.0.104.gtf')
gtf_df=as.data.frame(gtf)
annot = gtf_df
dim(annot)
names(annot)

# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$transcript_id)
probes2annot
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# 15520 transcripts do not have annotation using ensembl annotation. 
# This is a byproduct of the poor public transcriptome annotation
# for rats. Use the 2020 RTR transcriptome used to generate index transcriptome instead!

gtf <- rtracklayer::import('/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/RTR.gtf')
gtf_df=as.data.frame(gtf)
annot = gtf_df
dim(annot)
names(annot)

# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$transcript_id)
probes2annot
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
#0 transcripts are missing annotation

#=====================================================================================
#
#  Code chunk 23
#
#=====================================================================================


# Create the starting data frame
geneInfo0 = data.frame(transcript_id = probes,
                       gene_name = annot$gene_name[probes2annot],
                       gene_ID = annot$gene_id[probes2annot],
                       gene_name = annot$gene_name[probes2annot],
                       ref_gene_ID = annot$ref_gene_id[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for treatment
modOrder = order(-abs(cor(MEs, treatment, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.treatment));
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk 24
#
#=====================================================================================

#Filter down to significant modules:
geneInfo1 <- geneInfo[which(geneInfo$moduleColor=="turquoise" | geneInfo$moduleColor=="brown" | geneInfo$moduleColor=="thistle" | geneInfo$moduleColor=="purple" | geneInfo$moduleColor=="skyblue"),]
#7966 transcripts
write.csv(geneInfo1, file = "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_Male_PFF_sigmodules_geneInfo_062121.csv")

#Stratify dataset into significant modules:
geneInfo.turquoise <- geneInfo[which(geneInfo$moduleColor=="turquoise"),] #5337 transcripts
geneInfo.brown <- geneInfo[which(geneInfo$moduleColor=="brown"),] #1858 transcripts
geneInfo.thistle <- geneInfo[which(geneInfo$moduleColor=="thistle"),] #49 transcripts
geneInfo.purple <- geneInfo[which(geneInfo$moduleColor=="purple"),] #473 transcripts
geneInfo.skyblue <- geneInfo[which(geneInfo$moduleColor=="skyblue"),] #249 transcripts

#=====================================================================================
#
#  Code chunk 25
#
#=====================================================================================

####Finding hub genes

#Sort in descending order by MM.turquoise variable:
geneInfo.turquoise.sort <- geneInfo.turquoise[order(-abs(geneInfo.turquoise$MM.turquoise)),] #5337 genes
#Pull out top 10% by MM value (top 534 genes)
geneInfo.turquoise.sort.lim <- geneInfo.turquoise.sort[c(1:534),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.turquoise.sort.lim1 <- geneInfo.turquoise.sort.lim[c(1:16)]
turquoise_hub <- geneInfo.turquoise.sort.lim1
write.csv(turquoise_hub, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_turquoisemodule_hub_genes_Male_062121.csv")

#Sort in descending order by MM.brown variable:
geneInfo.brown.sort <- geneInfo.brown[order(-abs(geneInfo.brown$MM.brown)),] #1858 genes
#Pull out top 10% by MM value (top 186 genes)
geneInfo.brown.sort.lim <- geneInfo.brown.sort[c(1:186),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.brown.sort.lim1 <- geneInfo.brown.sort.lim[c(1:16)]
brown_hub <- geneInfo.brown.sort.lim1
write.csv(brown_hub, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_brownmodule_hub_genes_Male_062121.csv")

#Sort in descending order by MM.thistle variable:
geneInfo.thistle.sort <- geneInfo.thistle[order(-abs(geneInfo.thistle$MM.thistle)),] #49 genes
#Pull out top 10% by MM value (top 5 genes)
geneInfo.thistle.sort.lim <- geneInfo.thistle.sort[c(1:5),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.thistle.sort.lim1 <- geneInfo.thistle.sort.lim[c(1:16)]
thistle_hub <- geneInfo.thistle.sort.lim1
write.csv(thistle_hub, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_thistlemodule_hub_genes_Male_062121.csv")

#Sort in descending order by MM.purple variable:
geneInfo.purple.sort <- geneInfo.purple[order(-abs(geneInfo.purple$MM.purple)),] #473 genes
#Pull out top 10% by MM value (top 47 genes)
geneInfo.purple.sort.lim <- geneInfo.purple.sort[c(1:47),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.purple.sort.lim1 <- geneInfo.purple.sort.lim[c(1:16)]
purple_hub <- geneInfo.purple.sort.lim1
write.csv(purple_hub, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_purplemodule_hub_genes_Male_062121.csv")

#Sort in descending order by MM.skyblue variable:
geneInfo.skyblue.sort <- geneInfo.skyblue[order(-abs(geneInfo.skyblue$MM.skyblue)),] #249 genes
#Pull out top 10% by MM value (top 25 genes)
geneInfo.skyblue.sort.lim <- geneInfo.skyblue.sort[c(1:25),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.skyblue.sort.lim1 <- geneInfo.skyblue.sort.lim[c(1:16)]
skyblue_hub <- geneInfo.skyblue.sort.lim1
write.csv(skyblue_hub, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_skybluemodule_hub_genes_Male_062121.csv")


#Now compare lists of hub genes to the list of differentially expressed genes from sleuth:

#Read in list of male DGE results
male.dge <- read.csv("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/RNAseq_Sleuth_DGE_Male_TestTable_v3_060921.csv")
head(male.dge)
#Now cut down to only those with qvalue < 0.2
male.dge.q <- male.dge[which(male.dge$qval < 0.2),] #2534 observations
#Create transcript_id variable from "target_id" variable in female dge results:
male.dge.q$transcript_id <- male.dge.q$target_id
#Now merge female dge results from sleuth with the blue module:
turquoise_hub.dge <- merge(turquoise_hub, male.dge.q, by = "transcript_id")
#453 transcripts out of 534
brown_hub.dge <- merge(brown_hub, male.dge.q, by = "transcript_id")
#120 transcripts out of 186
thistle_hub.dge <- merge(thistle_hub, male.dge.q, by = "transcript_id")
#1 transcript out of 5
purple_hub.dge <- merge(purple_hub, male.dge.q, by = "transcript_id")
#6 transcripts out of 47
skyblue_hub.dge <- merge(skyblue_hub, male.dge.q, by = "transcript_id")
#0 transcripts out of 25

#Now join together the overlapping transcripts into a single dataframe:
module_hub1 <- rbind(turquoise_hub.dge, brown_hub.dge)
module_hub2 <- rbind(module_hub1, thistle_hub.dge)
module_hub3 <- rbind(module_hub2, purple_hub.dge)
all_hub.dge <- module_hub3

#Now write csv of all hub genes that are also differentially expressed (FDR<0.2 cutoff)
write.csv(all_hub.dge, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Male_062121.csv")
all_hub.dge <- read.csv("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Male_062121.csv")

## 6/3/21 - Next step: Try enrichr for GO analysis??
##Try performing enrichGO from clusterprofiler in separate R code.
##Try performing STRING on all genes.


## 6/21/21 - Compare list of hub+significant genes between males and females.

all_hub.dge.m <- read.csv("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Male_062121.csv")
#580 transcripts
all_hub.dge.f <- read.csv("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Female_051921.csv")
#331 transcripts

#Now try merging by shared ID:
all_hub.dge.merge <- merge(all_hub.dge.m, all_hub.dge.f, by = "transcript_id")
#22 transcripts show overlap between male and female WGCNA results

#Now write results as .csv:
write.csv(all_hub.dge.merge, "/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Male+Female_Overlap_062121.csv")

#Which of these overlapping transcripts show the same direction?

all_hub.dge.merge.up <- all_hub.dge.merge[which(all_hub.dge.merge$b.x > 0 & all_hub.dge.merge$b.y > 0),] 
#12 transcripts

all_hub.dge.merge.down <- all_hub.dge.merge[which(all_hub.dge.merge$b.x < 0 & all_hub.dge.merge$b.y < 0),] 
#5 transcripts



