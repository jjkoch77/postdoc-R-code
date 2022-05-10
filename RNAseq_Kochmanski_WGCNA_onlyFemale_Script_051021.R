#5/5/21 - Code for viewing existing sleuth results saved as .RData files

#First, you will need to load the R packages.
install.packages("remotes")
install.packages("WGCNA")
source("https://bioconductor.org/biocLite.R")
install.packages("BiocManager")
BiocManager::install(c("pachterlab/sleuth"))
BiocManager::install("WGCNA")
BiocManager::install("tidyverse")
library(WGCNA)
library(sleuth)
library(tidyverse)

# Next, load the previously saved RData files of RNA-seq sleuth results. Make sure to change /path/to/directory/ to correspond to the location of the
# files on your local computer.

#Log-likelihood ratio test results (differential expression hypothesis testing):
load("/Users/JoeKochmanski/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/AnnaStoll_DGE_sleuth_object_female_LRT_040921.RData")
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

colnames(so.F.tpm)
head(so.F.tpm)

#Note: We suggest removing features whose counts are consistently low (for example, removing all features that have a count of less than say 10 in more than 90% of the samples) 
#because such low-expressed features tend to reflect noise and correlations based on counts that are mostly zero aren't really meaningful. 
#The actual thresholds should be based on experimental design, sequencing depth and sample counts.

#Here, lets filter out all transcripts where at least one sample has less than 1 tpm:
so.F.tpm1 <- as.data.frame(so.F.tpm)
so.F.tpm2 <- so.F.tpm1[which(so.F.tpm1$sampF12P_quant > 1 & so.F.tpm1$sampF13P_quant > 1 & so.F.tpm1$sampF14P_quant > 1 & so.F.tpm1$sampF16P_quant > 1 & 
                                so.F.tpm1$sampF18C_quant > 1 & so.F.tpm1$sampF19C_quant > 1 & so.F.tpm1$sampF20C_quant > 1 & so.F.tpm1$sampF22C_quant > 1), ]
#27349 (out of 113100 total obs) transcripts have > 1 tpm for all samples.


#Now follow the WGCNA tutorial for consensus analysis of multiple, related datasets (e.g. male/female)
#Tutorial available here: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-DataInput.R

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
femData = so.F.tpm2;
femData = as.data.frame(femData)
# Take a quick look at what is in the data sets (caution, longish output):
dim(femData)
names(femData)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


datExpr0 = as.data.frame(t(femData));
names(datExpr0) = rownames(femData);
rownames(datExpr0) = names(femData);

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
abline(h = 60000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 60000, minSize = 4)
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


traitData = read.csv("/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/RNAseq_Sortwell_Meta_WGCNA_051021.csv");
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData
allTraits = allTraits[, c(1,2,5)];
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$sample);
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


save(datExpr, datTraits, file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/FemaleBrain-01-dataInput.RData")

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
lnames = load(file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/FemaleBrain-01-dataInput.RData");
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


#=====================================================================================
#
#  Code chunk 12
#
#=====================================================================================


net = blockwiseModules(datExpr, power = 10,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM", 
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
     file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/FemaleBrain-02-networkConstruction-auto.RData")


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
lnames = load(file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/FemaleBrain-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/FemaleBrain-02-networkConstruction-auto.RData");
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
#n=4 out of 76 modules with p-value <0.05:
#MEblue	p-value=5.756849e-05
#MEred	p-value=2.949849e-03
#MEyellow	p-value=1.193997e-02
#MEsalmon4	p-value=3.666340e-02

#The blue module is the most significant module! Then red, yellow, and salmon4, respectively.

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


module = "blue"
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

module = "red"
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

module = "yellow"
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

module = "salmon4"
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


names(datExpr)[moduleColors=="blue"]
names(datExpr)[moduleColors=="red"]
names(datExpr)[moduleColors=="yellow"]
names(datExpr)[moduleColors=="salmon4"]


#=====================================================================================
#
#  Code chunk 22
#
#=====================================================================================


#annot = read.csv(file = "GeneAnnotation.csv");
gtf <- rtracklayer::import('/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/Rattus_norvegicus.Rnor_6.0.104.gtf')
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
# 15437 transcripts do not have annotation using ensembl annotation. 
# This is a byproduct of the poor public transcriptome annotation
# for rats. Use the 2020 RTR transcriptome used to generate index transcriptome instead!

gtf <- rtracklayer::import('/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/RTR.gtf')
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
#0 transcripts do not have annotation

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
geneInfo1 <- geneInfo[which(geneInfo$moduleColor=="blue" | geneInfo$moduleColor=="red" | geneInfo$moduleColor=="yellow" | geneInfo$moduleColor=="salmon4"),]
#6419 transcripts
write.csv(geneInfo1, file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_Female_PFF_sigmodules_geneInfo_051221.csv")

#Stratify dataset into significant modules:
geneInfo.blue <- geneInfo[which(geneInfo$moduleColor=="blue"),] #3407 transcripts
geneInfo.red <- geneInfo[which(geneInfo$moduleColor=="red"),] #1185 transcripts
geneInfo.yellow <- geneInfo[which(geneInfo$moduleColor=="yellow"),] #1724 transcripts
geneInfo.salmon4 <- geneInfo[which(geneInfo$moduleColor=="salmon4"),] #103 transcripts

write.csv(geneInfo1, file = "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_Female_PFF_sigmodules_geneInfo_051221.csv")

#=====================================================================================
#
#  Code chunk 25
#
#=====================================================================================

####Finding hub genes

#Sort in descending order by MM.blue variable:
geneInfo.blue.sort <- geneInfo.blue[order(-abs(geneInfo.blue$MM.blue)),] #3407 genes
#Pull out top 10% by MM value (top 341 genes)
geneInfo.blue.sort.lim <- geneInfo.blue.sort[c(1:341),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.blue.sort.lim1 <- geneInfo.blue.sort.lim[c(1:16)]
blue_hub <- geneInfo.blue.sort.lim1
write.csv(blue_hub, "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_bluemodule_hub_genes_Female_051921.csv")

#Sort in descending order by MM.red variable:
geneInfo.red.sort <- geneInfo.red[order(-abs(geneInfo.red$MM.red)),] #1185 genes
#Pull out top 10% by MM value (top 119 genes)
geneInfo.red.sort.lim <- geneInfo.red.sort[c(1:119),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.red.sort.lim1 <- geneInfo.red.sort.lim[c(1:16)]
red_hub <- geneInfo.red.sort.lim1
write.csv(red_hub, "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_redmodule_hub_genes_Female_051921.csv")

#Sort in descending order by MM.yellow variable:
geneInfo.yellow.sort <- geneInfo.yellow[order(-abs(geneInfo.yellow$MM.yellow)),] #1724 genes
#Pull out top 10% by MM value (top 172 genes)
geneInfo.yellow.sort.lim <- geneInfo.yellow.sort[c(1:172),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.yellow.sort.lim1 <- geneInfo.yellow.sort.lim[c(1:16)]
yellow_hub <- geneInfo.yellow.sort.lim1
write.csv(yellow_hub, "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_yellowmodule_hub_genes_Female_051921.csv")

#Sort in descending order by MM.salmon4 variable:
geneInfo.salmon4.sort <- geneInfo.salmon4[order(-abs(geneInfo.salmon4$MM.salmon4)),] #103 genes
#Pull out top 10% by MM value (top 10 genes)
geneInfo.salmon4.sort.lim <- geneInfo.salmon4.sort[c(1:10),]
#Now limit variables to only relevant columns for the sig modules:
geneInfo.salmon4.sort.lim1 <- geneInfo.salmon4.sort.lim[c(1:16)]
salmon4_hub <- geneInfo.salmon4.sort.lim1
write.csv(salmon4_hub, "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_salmon4module_hub_genes_Female_051921.csv")

#Now compare lists of hub genes to the list of differentially expressed genes from sleuth:

#Read in list of female DGE results
female.dge <- read.csv("/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/RNAseq_Sleuth_DGE_Female_TestTable_v2_040921.csv")
head(female.dge)
#Now cut down to only those with qvalue < 0.2
female.dge.q <- female.dge[which(female.dge$qval < 0.2),] #1419 observations
#Create transcript_id variable from "target_id" variable in female dge results:
female.dge.q$transcript_id <- female.dge.q$target_id
#Now merge female dge results from sleuth with the blue module:
blue_hub.dge <- merge(blue_hub, female.dge.q, by = "transcript_id")
#248 transcripts
red_hub.dge <- merge(red_hub, female.dge.q, by = "transcript_id")
#63 transcripts
yellow_hub.dge <- merge(yellow_hub, female.dge.q, by = "transcript_id")
#20 transcripts
salmon4_hub.dge <- merge(salmon4_hub, female.dge.q, by = "transcript_id")
#0 transcripts

#Now join together the overlapping transcripts into a single dataframe:
module_hub1 <- rbind(blue_hub.dge, red_hub.dge)
module_hub2 <- rbind(module_hub1, yellow_hub.dge)
module_hub3 <- rbind(module_hub2, salmon4_hub.dge)
all_hub.dge <- module_hub3

#Now write csv of all hub genes that are also differentially expressed (FDR<0.2 cutoff)
write.csv(all_hub.dge, "/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Female_051921.csv")
all_hub.dge <- read.csv("/Users/jjkoch/Dropbox/Bernstein_Lab/Random TSMM Stuff/AnnaStoll_RNAseq/WGCNA_allsighub_siggenes_Female_051921.csv")

## 5/12/21 - Next step: Try enrichr for GO analysis??
##5/19/21 - Performed enrichGO from clusterprofiler in separate R code.

