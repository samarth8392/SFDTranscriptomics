###########################################################################
###                           Samarth Mathur, PhD                       ###
###                         The Ohio State University                   ###
###                                                                     ###
###     Date Created: 04/19/21                  Last Modified: 06/23/21 ###
###########################################################################
###########################################################################
###                   WGCNA_kidney.R       		                        ###
###########################################################################

### Weighted gene co-expression network analysis (WGCNA) for RNA-seq data ###
### Largely based on WGCNA tutorials at: 
### https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

#### PREREQUISITES #####

# Install/load libraries

#install.packages("BiocManager", lib="/scratch/bell/mathur20/osu/rnaseq/Rlibs", repos="http://cran.us.r-project.org")
#BiocManager::install("DESeq2", lib="/scratch/bell/mathur20/osu/rnaseq/Rlibs")
#BiocManager::install("WGCNA", lib="/scratch/bell/mathur20/osu/rnaseq/Rlibs")
#install.packages("dynamicTreeCut", lib="/scratch/bell/mathur20/osu/rnaseq/Rlibs", repos="http://cran.us.r-project.org")
#install.packages("fastcluster", lib="/scratch/bell/mathur20/osu/rnaseq/Rlibs", repos="http://cran.us.r-project.org")


library("DESeq2", lib.loc ="/scratch/bell/mathur20/osu/rnaseq/Rlibs")
library("dynamicTreeCut", lib.loc ="/scratch/bell/mathur20/osu/rnaseq/Rlibs")
library("fastcluster", lib.loc ="/scratch/bell/mathur20/osu/rnaseq/Rlibs")
library("WGCNA", lib.loc ="/scratch/bell/mathur20/osu/rnaseq/Rlibs")

# 

setwd("/scratch/bell/mathur20/osu/rnaseq/WGCNA/kidney")
load ("/scratch/bell/mathur20/osu/rnaseq/DESeq2/Final/Kidney/through_logxform_kidney.RData")

#### STEP1: Data input and cleanup #####
new.dds_kidney <- estimateSizeFactors(dds_kidney)
counts_kidney <- counts(new.dds_kidney, normalized=TRUE)

options(stringsAsFactors = FALSE)
dim(counts_kidney)  # 107167 transcripts, 11 samples

head(rownames(counts_kidney))

## check for genes and samples with too many missing values
gsg <- goodSamplesGenes(counts_kidney, verbose = 3)
gsg$allOK ## <-- this should return "TRUE" ==> if it doesn't, run below loop

# if (!gsg$allOK)
# {
#   ## this is where genes that weren't "good" were filtered out.
#   ## optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(counts)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(counts)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   counts = counts[gsg$goodSamples, gsg$goodGenes]
# }

## how many genes have at least 8 samples with non-zero counts?
#dim(counts_kidney[which(rowSums(counts_kidney) != 0),]) # 51301

dim(counts_kidney[rowSums(counts_kidney==0) <= 8,]) # 47522
counts_kidney <- counts_kidney[rowSums(counts_kidney==0) <= 8,]

sampleTree = hclust(dist(t(counts_kidney)), method = "average")

pdf(file = "sampleClustering.pdf", width = 12, height = 9)
par(cex = 1.5)
par(mar = c(2,4.5,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
      cex.axis = 1, cex.main = 1)
#abline(h=3.5e5, col="red")
dev.off()  # No outliers

#clust = cutreeStatic(sampleTree, cutHeight = 3.5e5, minSize = 9) # All remaining 9 out of 10 samples are in cluster 1
#table(clust) 
#keepSamples = (clust==1)
#datExpr = counts_kidney[,keepSamples]
datExpr <- t(counts_kidney)
nGenes = ncol(datExpr) # 47522
nSamples = nrow(datExpr) # 11

# The variable datExpr now contains the expression data ready for network analysis.

collectGarbage()

##### Step2: Network Construction #####

enableWGCNAThreads(n=64)

# choose a set of soft-thresholding powers to start with
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
## call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
## plot the results

## scale-free topology fit index as a function of the soft-thresholding power

pdf(file = "threshold.pdf", width = 9, height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab=expression(paste("Scale Free Topology Model Fit,signed ",R^2)),type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
## mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()  # We choose power=8


# Co-expression similarity and adjacency
# All genes at once (no block usage)

softPower = 8
adjacency_kidney = adjacency(datExpr, power = softPower, type = "signed")

# Topological Overlap Matrix (TOM)
TOM_kidney = TOMsimilarity(adjacency_kidney, TOMType = "signed", verbose = 5)
dissTOM_kidney= 1-TOM_kidney

# Clustering using TOM

# Call the hierarchical clustering function
geneTree_kidney = hclust(as.dist(dissTOM_kidney), method = "average")

# Plot the resulting clustering tree (dendrogram)
pdf(file = "geneTree_Kidney.pdf", width = 12, height = 9)
plot(geneTree_kidney, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity (Kidney)",
labels = FALSE, hang = 0.04)
dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree_kidney, distM = dissTOM_kidney,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
table(dynamicMods) # 179 modules

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "geneTree_Kidney_dynamicMod.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree_kidney, dynamicColors, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors (Kidney)")
dev.off()

#save.image("through_Step2_NetworkCon_Kidney.RData")
#  Merging of modules whose expression profiles are very similar

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
pdf(file = "eigengene_cluster_kidney.pdf", width = 30, height = 9)
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")
#abline(h=MEDissThres, col = "red")
dev.off()

# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
# Trying different cut heights
MEDissThres = c(0.1,0.25,0.4,0.75)

# Call an automatic merging function
merge1 = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres[1], verbose = 3)
merge2 = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres[2], verbose = 3)
merge3 = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres[3], verbose = 3)
merge4 = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres[4], verbose = 3)

# The merged module colors
mergedColors1 = merge1$colors
mergedColors2 = merge2$colors
mergedColors3 = merge3$colors
mergedColors4 = merge4$colors
# Eigengenes of the new merged modules:
mergedMEs1 = merge1$newMEs
mergedMEs2 = merge2$newMEs
mergedMEs3 = merge3$newMEs
mergedMEs4 = merge4$newMEs

pdf(file = "geneTree_Kidney_MergedMod.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree_kidney, cbind(dynamicColors, mergedColors1,mergedColors2,mergedColors3,mergedColors4),
c("Dynamic Tree Cut", "Merged dynamic 0.1", "Merged dynamic 0.25", "Merged dynamic 0.4", "Merged dynamic 0.75"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off() # Best is 0.25

pdf(file = "geneTree_Kidney_MergedMod.pdf", width = 12, height = 6)
plotDendroAndColors(geneTree_kidney, cbind(mergedColors2),
c("Modules (N=62)"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off() # Best is 0.25


table(mergedColors2)
length(table(mergedColors2)) # N=62

save.image("through_Step2_NetworkCon_Kidney.RData")

##### Step3: Relating modules to external traits #####

load("through_Step2_NetworkCon_Kidney.RData")

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Load experimental design (trait data)
# Male = 1, Female = 2; 20C = 20 , 26C = 26; Control = 0, Treated = 1 
all_traits_kidney <- read.csv("Traits_Kidney.csv")
#all_traits_kidney <- all_traits_kidney[c(2:10),] # Remove outlier
all_traits_kidney <- all_traits_kidney[,c(1,2,4,5)] # Remove sex

samples = rownames(datExpr)
traitRows = match(samples, all_traits_kidney$Sample_ID)
datTraits = all_traits_kidney[traitRows, -c(1,2)]
rownames(datTraits) = all_traits_kidney[traitRows, 1]

# Rename to moduleColors
moduleColors = mergedColors2
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(53))
moduleLabels = match(mergedColors2, colorOrder)-1
MEs = mergedMEs2


# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr,moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#pdf(file = "Module-trait_relationship_kidney.pdf", width = 9, height = 12)
pdf(file = "Module-trait_relationship_kidneyA.pdf", width = 20, height = 5)
#par(mar = c(6, 8.5, 3, 3))
par(mar = c(12, 6, 6, 8.5))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = t(moduleTraitCor),
yLabels = names(datTraits),
xLabels = substring(names(MEs), 3),
xSymbols = substring(names(MEs), 3),
colorLabels = FALSE,
colors = blueWhiteRed(53, gamma = 1, endSaturation = 1),
textMatrix = t(textMatrix),
setStdMargins = FALSE,
cex.text = 0.5,
zlim = c(-1,1),
main = paste("Module-trait relationships"))
dev.off()

# Gene relationship to trait and important modules
# Gene Significance and Module Membership

# Define variable treatment containing the treatment column of datTrait
treatment = as.data.frame(datTraits$Treatment)
names(treatment) = "Treatment"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, treatment, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(treatment), sep="")
names(GSPvalue) = paste("p.GS.", names(treatment), sep="")

# Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
# membership in interesting modules. 
# we look at the  module that has the highest association with treatment
# antiquewhite1 (p=0.001) darkgreen (p=0.02) burlywood (p=0.03) slateblue1 (p=0.03) magenta2 (p=0.04) 
# We plot a scatterplot of Gene Significance vs. Module Membership each module

### Step 4: Visualizing the gene network ###

# Only for 400 random genes
nSelect = 400
# For reproducibility, we set the random seed
set.seed(44)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM_kidney[select, select]
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window
pdf(file = "GeneNetwork_heatmap_kidney.pdf", width = 9, height = 9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^10
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, terrainColors=T,
	main = "Network heatmap plot, (4,00 random genes)")
dev.off()
#save.image("through_Step2_NetworkCon_Kidney.RData")

# Visualizing the network of eigengenes

# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate treatment from the clinical traits
treatment = as.data.frame(datTraits$Treatment)
temp = as.data.frame(datTraits$Temp)
names(treatment) = "treatment"
names(temp) = "temp"
# Add the treatment to existing module eigengenes
MET = orderMEs(cbind(MEs, treatment))
# Plot the relationships among the eigengenes and the trait
pdf(file = "GeneNetwork_eigenTreat_kidney.pdf", width = 9, height = 6)
par(cex = 0.9)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
plotHeatmaps = FALSE)
dev.off()
# Plot the dendrogram
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf(file = "GeneNetwork_eigenheatTreat_kidney.pdf", width = 9, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(6,8,5,7),
plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

### Step 5: Getting list of gene IDs ###
sig.modules <- moduleTraitPvalue[which(moduleTraitPvalue[,3] < 0.045),]
modulelist <- rownames(sig.modules)
modulelist <- substring(modulelist, 3)
# adding more modules based on adajency matrix (i.e. closeness of modules)

moremodules <- c("blue4","indianred4","paleturquoise4","brown3","cornflowerblue")

modulelist <- c(modulelist,moremodules)

pdf(file = "GSvMM_Kidney.pdf", width = 12, height = 9)
par(mfrow = c(3,4))
for (i in 1:12)
{
	column = match(modulelist[i], modNames)
	moduleGenes = moduleColors==modulelist[i]
	verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
	abs(geneTraitSignificance[moduleGenes, 1]),
	xlab = paste("Module Membership in", modulelist[i], "module"),
	ylab = "Gene significance for treatment",
	main = paste("Module membership vs. gene significance\n"),
	cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col="black", bg = modulelist[i], pch=21)
}
dev.off()

# Summary output of network analysis results

finalmodules <- modulelist[c(1,2,3,4,6,7,8,9,10,11)]

moduleGenes1 <- datExpr[,moduleColors==finalmodules[1]]
moduleGenes2 <- datExpr[,moduleColors==finalmodules[2]]
moduleGenes3 <- datExpr[,moduleColors==finalmodules[3]]
moduleGenes4 <- datExpr[,moduleColors==finalmodules[4]]
moduleGenes5 <- datExpr[,moduleColors==finalmodules[5]]
moduleGenes6 <- datExpr[,moduleColors==finalmodules[6]]
moduleGenes7 <- datExpr[,moduleColors==finalmodules[7]]
moduleGenes8 <- datExpr[,moduleColors==finalmodules[8]]
moduleGenes9 <- datExpr[,moduleColors==finalmodules[9]]
moduleGenes10 <- datExpr[,moduleColors==finalmodules[10]]

geneinModules <- c(dim(moduleGenes1)[2],dim(moduleGenes2)[2],dim(moduleGenes3)[2],
	dim(moduleGenes4)[2],dim(moduleGenes5)[2],dim(moduleGenes6)[2], dim(moduleGenes7)[2],
	dim(moduleGenes8)[2],dim(moduleGenes9)[2],dim(moduleGenes10)[2])
geneinModules
# [1] 2684  145   96   97  123  138 1579  287   49  152

sum(geneinModules) # N=5350
list.sig.modules <- finalmodules

module.gene.IDs <- data.frame()
i <- 1
for (mod in finalmodules) {
  tmp.gene.IDs <- data.frame()
  mod.genes <- dimnames(datExpr)[[2]][moduleColors==mod]
  tmp.gene.IDs[i:(length(mod.genes)),1] <- mod
  tmp.gene.IDs[i:(length(mod.genes)),2] <- mod.genes
  module.gene.IDs <- rbind(module.gene.IDs, tmp.gene.IDs)
}

sum(table(module.gene.IDs)) # 5350
write.csv(module.gene.IDs, "Kidney_Module_gene_IDs.csv")


# List of modules significantly associated with Treatment, temperature, or sex #
head(moduleTraitCor)
head(moduleTraitPvalue)

all(rownames(moduleTraitCor) == rownames(moduleTraitPvalue))
all(colnames(moduleTraitCor) == colnames(moduleTraitPvalue))

p.vals <- data.frame(module=character(),
                     n.genes=numeric(),
                     pred.var=character(),
                     cor=numeric(),
                     p.value=numeric(),
                     stringsAsFactors=FALSE)
temp <- data.frame(module=character(),
                   n.genes=numeric(),
                   pred.var=character(),
                   cor=numeric(),
                   p.value=numeric(),
                   stringsAsFactors=FALSE)

for(row in 1:nrow(moduleTraitPvalue)) {
  for (col in 1:ncol(moduleTraitPvalue)) {
    if (moduleTraitPvalue[row,col]<0.045) {
      temp.mod <- substring(rownames(moduleTraitPvalue)[row],3)
      temp[1,1] <- rownames(moduleTraitPvalue)[row]
      temp[1,3] <- colnames(moduleTraitPvalue)[col]
      temp[1,4] <- moduleTraitCor[row,col]
      temp[1,5] <- moduleTraitPvalue[row,col]
      temp[1,2] <- length(colnames(datExpr)[moduleColors==temp.mod])
      p.vals <- rbind(p.vals, temp)
    }
  }
}

write.csv(p.vals,"Kidney_Sig_modules_cors_pvals.csv", row.names=F)

# Retain genes with geneIDs

modulegenes <- read.csv("Kidney_Module_gene_IDs.csv", header=T)
modulegenes <- modulegenes[-which(is.na(modulegenes$GeneID)),]

dim(modulegenes) #  2643    4
table(modulegenes$Module)
#         blue4         brown3         green2     indianred4        magenta 
#           921             90             71            173           1129 
#      magenta4    navajowhite paleturquoise4        salmon1       thistle4 
#            78             61             20             45             55 

write.csv(modulegenes,"Kidney_Module_gene_IDs_noNA.csv", row.names=F)

save.image("through_Step2_NetworkCon_Kidney.RData")


### Step6: Visualize gene networks ####

load("through_Step2_NetworkCon_Kidney.RData")

# Exporting to VisANT

finalmodules

# Read in the annotation file
annot = read.csv(file = "/scratch/bell/mathur20/osu/rnaseq/ref/annotation/CroVir_GeneAnnotation.csv")




# Select module
module = finalmodules[4]
# Select module probes
probes = colnames(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM_kidney[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into an edge list file VisANT can read

#vis = exportNetworkToVisANT(modTOM,
#file = paste("VisANTInput-", module, ".txt", sep=""),
#weighted = TRUE,
#threshold = 0,
#probeToGene = data.frame(annot$GeneID, annot$Homolog_Anolis) )

# we can restrict the genes in the output to say the 30 top hub genes in the
# module
#nTop = 30
#MConn = softConnectivity(datExpr[, modProbes])
#top = (rank(-IMConn) <= nTop)
#vis = exportNetworkToVisANT(modTOM[top, top],
#file = paste("VisANTInput-", module, "-top30_Kidney.txt", sep=""),
#weighted = TRUE,
#threshold = 0)
#probeToGene = data.frame(annot$GeneID, annot$Homolog_Anolis) )

# Cytoscape
setwd("networks")
for (i in 1:length(finalmodules)) 
{
  modules = c(finalmodules[i])
  # Select module probes
  probes = colnames(datExpr)
  inModule = is.finite(match(moduleColors, modules));
  modProbes = probes[inModule];
  #modGenes = annot$Homolog_Anolis[match(modProbes, annot$substanceBXH)];
  # Select the corresponding Topological Overlap
  modTOM = TOM_kidney[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)

  cyt = exportNetworkToCytoscape(modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), "_kidney.txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), "_kidney.txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]);
}

# Clean up and modify networks


