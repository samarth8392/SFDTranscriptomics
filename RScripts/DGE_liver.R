###########################################################################
###                           Samarth Mathur, PhD                       ###
###                         The Ohio State University                   ###
###                                                                     ###
###     Date Created: 06/01/21                  Last Modified: 08/09/21 ###
###########################################################################
###########################################################################
###                   DGE_final_liver.R   		                          ###
###########################################################################

## Analysis based on : http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#### PREREQUISITES #####
library("apeglm")
library("ashr")
library("DESeq2")
library("ggplot2")
library("genefilter")
library("dendextend")
library("rtracklayer")
library("gplots")
library("ggrepel")
library("vsn")
library("hexbin")
library(dplyr)
library(VennDiagram)
library("pheatmap")
setwd("~/Documents/Postdoc_research/rnaseq/results/")

# Analyzing each kind of tissue separately (10 liver, 11 kidney, 10 skin)

####### ONLY LIVER (NO CV16) #########
#### Data loading and prep ####
gene.counts <- read.table("featurecounts/featureCount_counts_liver_controlfirst_no16.txt", sep="\t", header=T, row.names=1, check.names=F)
drop.vars <- names(gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
gene.counts <- gene.counts[!drop.vars]
head(gene.counts)

assembly <- readGFF("stringtie_all_merged.annotated.gtf")  

## read in CSV with sample information
col.data <- read.csv("Experimental_design_Liver_no16.csv", row.names = 1)

#Change the colnames in gene.counts to match experimental design
colnames(gene.counts) <- rownames(col.data)

## set temperature as a factor
col.data[,3] <- as.factor(col.data[,3])

## check to make sure that sample names are in the same order in the gene count table and the
## sample info table
all(rownames(col.data) %in% colnames(gene.counts))
all(rownames(col.data) == colnames(gene.counts))

#### DGE analysis #### 
# Construction of DESeqDataSet object and DE analysis
dds_liver <- DESeqDataSetFromMatrix(countData = gene.counts,
                              colData = col.data,
                              design = ~ Temp + Treatment + Temp:Treatment)

## renames genes in a DESeqDataSet to match gene names in the reference annotation 
## (i.e., replaces MSTRG wherever gene name is known)

gene_idx <- match(dds_liver@rowRanges@partitioning@NAMES, assembly$transcript_id)

## create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different 
## isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)

gene.name <- assembly$gene_name[gene_idx]
transcript.id <- assembly$transcript_id[gene_idx]
xloc <- assembly$xloc[gene_idx]
gene_names <- cbind(gene.name, transcript.id, xloc)

## adds unique gene names to dds
dds_liver@rowRanges@partitioning@NAMES <- paste(gene_names[,1],gene_names[,2],gene_names[,3], sep=" ")

## check to make sure that no gene names in dds are duplicates
which(duplicated(dds_liver@rowRanges@partitioning@NAMES))

# Run the statistical analysis
dds_liver <- DESeq(dds_liver)

## examine distribution of dispersion values
#plotDispEsts(dds)

## examine MA plot (normalized counts vs. log2fold changes)
#plotMA(res)

#### Data visualization #### 

## rlog() = transforms count data to the log2 scale in a way that minimizes differences between
## samples for genes with small counts, and which normalizes with respect to library size
rld_liver <- rlog(dds_liver)
#meanSdPlot(assay(rld_liver))
#head(assay(rld_liver))

save.image("Final/Liver/through_logxform_liver.RData")

# Heatmap of the sample-to-sample distances
library("RColorBrewer")
sampleDists <- dist(t(assay(rld_liver)), method = "euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld_liver$Individual, rld_liver$Temp, rld_liver$Treatment, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(file="Final/Liver/SampleCluster_Liver.pdf", width=8, height=9)
plot(hclust(sampleDists), cex=1.5, lwd=1.5)
#pheatmap(sampleDistMatrix,
#         clustering_distance_rows=sampleDists,
#         clustering_distance_cols=sampleDists,
#         col=colors)
dev.off()

# PCA to examine effects of treatment + temp
plot.treat.data <- plotPCA(rld_liver, intgroup = c("Treatment", "Temp"), returnData=TRUE, n=107167) ## uses n most variable genes -- not necessarily sig DEGs!
percentVar <- round(100*attr(plot.treat.data, "percentVar"))
my.colors <-c("#00BFC4","#F8766D")
pdf(file="Final/Liver/PCA_allgenes_Liver.pdf", width=9, height=6)
ggplot(plot.treat.data, aes(PC1, PC2, color=Treatment,shape=Temp)) +
  geom_point(size=4) +
  geom_text_repel(aes(label =name),
                  size = 4, max.overlaps = 32)+
  scale_color_manual(values=c(my.colors)) +
  theme_linedraw(base_size = 18) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
dev.off()

#### DEG Identification #### 
# Analysis based on : https://github.com/tavareshugo/tutorial_DESeq2_contrasts

#get the model matrix
mod_mat <- model.matrix(design(dds_liver), colData(dds_liver))
# Coefficients
resultsNames(dds_liver) # resultsNames(dds_liver)
#[1] "Intercept"                    "Temp_26C_vs_20C"              "Treatment_Treated_vs_Control" "Temp26C.TreatmentTreated"

# Output DEGs by co-efficients
temp <- results(dds_liver, contrast = list("Temp_26C_vs_20C"),
                           alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Con vs Treat (at 20C in liver)
treat <- results(dds_liver, contrast = list("Treatment_Treated_vs_Control"),
                 alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Con vs Treat (at 20C in liver)
inter <- results(dds_liver, contrast = list("Temp26C.TreatmentTreated"),
                 alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Con vs Treat (at 20C in liver)

# get shrunken log fold changes
temp_shrink <- lfcShrink(dds_liver,contrast = list("Temp_26C_vs_20C"), type='ashr', res = temp)
treat_shrink <- lfcShrink(dds_liver,contrast = list("Treatment_Treated_vs_Control"), type='ashr', res = treat)
inter_shrink <- lfcShrink(dds_liver,contrast = list("Temp26C.TreatmentTreated"), type='ashr', res = inter)

pdf(file="Final/Liver/Shrink_Liver_Coeff.pdf", width=9, height=6)
par(mfrow=c(1,3))
plot(
  x=temp$log2FoldChange,
  y=temp_shrink$log2FoldChange,pch=20,
  cex=.5,
  col=1+(temp_shrink$padj < 0.05),
  xlab="raw log2 fold change",
  ylab=expression(paste("shrunken log2 fold change (",beta[1], ")")),
  xlim=c(-15,15),
  ylim=c(-15,15),
  main="DGE due to Temperature (26C vs 20C)"
)
abline(0,1)

plot(
  x=treat$log2FoldChange,
  y=treat_shrink$log2FoldChange,pch=20,
  cex=.5,
  col=1+(treat_shrink$padj < 0.05),
  xlab="raw log2 fold change",
  ylab=expression(paste("shrunken log2 fold change (",beta[1], ")")),
  xlim=c(-15,15),
  ylim=c(-15,15),
  main="DGE due to Treatment (SFD vs Con)"
)
abline(0,1)

plot(
  x=inter$log2FoldChange,
  y=inter_shrink$log2FoldChange,pch=20,
  cex=.5,
  col=1+(inter_shrink$padj < 0.05),
  xlab="raw log2 fold change",
  ylab=expression(paste("shrunken log2 fold change (",beta[1], ")")),
  xlim=c(-15,15),
  ylim=c(-15,15),
  main="DGE due to Interaction (SFD and 26C)"
)
abline(0,1)
dev.off()

# No. of DGEs
sum(temp_shrink$padj < 0.05, na.rm=T) # N = 413
sum(treat_shrink$padj < 0.05, na.rm=T) # N = 1175
sum(inter_shrink$padj < 0.05, na.rm=T) # N = 236

#Output results
resSig <- temp_shrink[which(temp_shrink$padj < 0.05),]
write.csv(as.data.frame(resSig), file="Final/Liver/DGE_Temp_Liver.csv", quote = F, row.names = T)
resSig <- treat_shrink[which(treat_shrink$padj < 0.05),]
write.csv(as.data.frame(resSig), file="Final/Liver/DGE_Treat_Liver.csv", quote = F, row.names = T)
resSig <- inter_shrink[which(inter_shrink$padj < 0.05),]
write.csv(as.data.frame(resSig), file="Final/Liver/DGE_Inter_Liver.csv", quote = F, row.names = T)


#Venn diagrams
res1 <- read.csv("Final/Liver/DGE_Temp_Liver.csv", header = T)
res2 <- read.csv("Final/Liver/DGE_Treat_Liver.csv", header = T)
res3 <- read.csv("Final/Liver/DGE_Inter_Liver.csv", header = T)

venn.diagram(list("Temp (26C vs 20C)"=res1$X, "Treat (SFD vs Con)"=res2$X, "Temp x Treat"=res3$X), fill = c("yellow","cyan","green"), 
             height = 4500, width = 4500, alpha=c(0.5,0.5,0.5),
             cex = 1, filename="Final/Liver/Venn_Liver_Coeff.png")


# # Output DEGs by contrast
# define coefficient vectors for each contrast
con_20 <- colMeans(mod_mat[dds_liver$Treatment == "Control" & dds_liver$Temp == "20C", ])
tre_20 <- colMeans(mod_mat[dds_liver$Treatment == "Treated" & dds_liver$Temp == "20C", ])
con_26 <- colMeans(mod_mat[dds_liver$Treatment == "Control" & dds_liver$Temp == "26C", ])
tre_26 <- colMeans(mod_mat[dds_liver$Treatment == "Treated" & dds_liver$Temp == "26C", ])

res1_con_tre_20 <- results(dds_liver, contrast = con_20 - tre_20, 
                           alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Con vs Treat (at 20C in liver)
res2_con_tre_26 <- results(dds_liver, contrast = con_26 - tre_26, 
                           alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Con v Treat (at 26C in liver)
res3_inter <- results(dds_liver, contrast = (con_20 - con_26) - (tre_20 - tre_26)) # Temp x Treat interaction (Con vs Treat (at 26C in liver))
# Interaction between Treatment and Temperature (i.e. do Con and Treated respond differently to the higher Temperature):

#res4_con_20_26 <- results(dds_liver, contrast = con_26 - con_20, 
#                          alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Con at lower temp (at 20C in liver)

res5_tre_20_26 <- results(dds_liver, contrast = (con_26 - con_20)-(tre_26 - tre_20), 
                          alpha=0.05, pAdjustMethod = "fdr", independentFiltering = TRUE) # Extra DEGs between treated due to temp


# get shrunken log fold changes
res1_shrink_b1 <- lfcShrink(dds_liver,contrast = con_20 - tre_20, type='ashr', res = res1_con_tre_20)
res2_shrink_b2 <- lfcShrink(dds_liver,contrast = con_26 - tre_26, type='ashr', res=res2_con_tre_26)
res3_shrink_b3 <- lfcShrink(dds_liver,contrast = (con_20 - con_26) - (tre_20 - tre_26), type='ashr', res=res3_inter)
#res4_shrink_b4 <- lfcShrink(dds_liver,contrast = con_26 - con_20, type='ashr', res=res4_con_20_26)
res5_shrink_b5 <- lfcShrink(dds_liver,contrast = (con_26 - con_20)-(tre_26 - tre_20), type='ashr', res=res5_tre_20_26)

# plot the shrunken log2 fold changes against the raw changes:
pdf(file="Final/Liver/Shrink_Liver_Contrast2.pdf", width=9, height=6)
par(mfrow=c(1,2))
plot(
  x=res3_inter$log2FoldChange,
  y=res3_shrink_b3$log2FoldChange,pch=20,
  cex=.5,
  col=1+(res3_shrink_b3$padj < 0.05),
  xlab="raw log2 fold change",
  ylab=expression(paste("shrunken log2 fold change")),
  xlim=c(-15,15),
  ylim=c(-15,15),
  main="DGE due to Interaction (Temp X Treat)"
)
abline(0,1)

plot(
  x=res5_tre_20_26$log2FoldChange,
  y=res5_shrink_b5$log2FoldChange,pch=20,
  cex=.5,
  col=1+(res5_shrink_b5$padj < 0.05),
  xlab="raw log2 fold change",
  ylab=expression(paste("shrunken log2 fold change")),
  xlim=c(-15,15),
  ylim=c(-15,15),
  main="DGE 26C vs 20C (for treated)"
)
abline(0,1)
dev.off()

# No. of DGEs
sum(res1_shrink_b1$padj < 0.05, na.rm=T) # N = 1175
sum(res2_shrink_b2$padj < 0.05, na.rm=T) # N = 132
sum(res3_shrink_b3$padj < 0.05, na.rm=T) # N = 231
#sum(res4_shrink_b4$padj < 0.05, na.rm=T) # N = 413
sum(res5_shrink_b5$padj < 0.05, na.rm=T) # N = 236

#Output results
resSig1 <- res1_shrink_b1[which(res1_shrink_b1$padj < 0.05),]
write.csv(as.data.frame(resSig1), file="Final/Liver/DGE_Con1_Liver.csv", quote = F, row.names = T)
resSig2 <- res2_shrink_b2[which(res2_shrink_b2$padj < 0.05),]
write.csv(as.data.frame(resSig2), file="Final/Liver/DGE_Con2_Liver.csv", quote = F, row.names = T)
resSig3 <- res3_shrink_b3[which(res3_shrink_b3$padj < 0.05),]
write.csv(as.data.frame(resSig3), file="Final/Liver/DGE_Con3_Liver.csv", quote = F, row.names = T)
resSig4 <- res4_shrink_b4[which(res4_shrink_b4$padj < 0.05),]
write.csv(as.data.frame(resSig4), file="Final/Liver/DGE_Con4_Liver.csv", quote = F, row.names = T)
resSig5 <- res5_shrink_b5[which(res5_shrink_b5$padj < 0.05),]
write.csv(as.data.frame(resSig5), file="Final/Liver/DGE_Con5_Liver.csv", quote = F, row.names = T)

#Venn diagrams
res1 <- read.csv("Final/Liver/DGE_Con1_Liver.csv", header = T)
res2 <- read.csv("Final/Liver/DGE_Con2_Liver.csv", header = T)
res3 <- read.csv("Final/Liver/DGE_Con3_Liver.csv", header = T)

venn.diagram(list("Con Vs Treated (20C)"=res1$X, "Con Vs Treated (26C)"=res2$X, "26C x Treat"=res3$X), fill = c("yellow","cyan","green"), 
             height = 4500, width = 4500, alpha=c(0.5,0.5,0.5),
             cex = 1, filename="Final/Liver/Venn_Liver_Constrast.png")

# PCA to examine effects of treatment + temp (DGEs)
p.cutoff <- 0.05
fc.cutoff <- 0

res1_topSigGenes <- temp_shrink[which(temp_shrink$padj < p.cutoff & abs(temp_shrink$log2FoldChange) >= fc.cutoff),]
res2_topSigGenes <- treat_shrink[which(treat_shrink$padj < p.cutoff & abs(treat_shrink$log2FoldChange) >= fc.cutoff),]
res3_topSigGenes <- inter_shrink[which(inter_shrink$padj < p.cutoff & abs(inter_shrink$log2FoldChange) >= fc.cutoff),]

num.genes1 <- res1_topSigGenes@nrows
length(rownames(res1_topSigGenes))
rld_liver.Sig1 <- rld_liver[rownames(rld_liver) %in% rownames(res1_topSigGenes)]

num.genes2 <- res2_topSigGenes@nrows
length(rownames(res2_topSigGenes))
rld_liver.Sig2 <- rld_liver[rownames(rld_liver) %in% rownames(res2_topSigGenes)]

num.genes3 <- res3_topSigGenes@nrows
length(rownames(res3_topSigGenes))
rld_liver.Sig3 <- rld_liver[rownames(rld_liver) %in% rownames(res3_topSigGenes)]

# heatmap of count matrix (DGE)

df <- as.data.frame(colData(dds_liver)[,c("Temp","Treatment")])
pdf(file="Final/Liver/CountMatrix_Liver_Temp.pdf", width=6, height=9)
pheatmap(assay(rld_liver.Sig1), cluster_rows=T, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df,
         clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()
pdf(file="Final/Liver/CountMatrix_Liver_Treat.pdf", width=6, height=9)
pheatmap(assay(rld_liver.Sig2), cluster_rows=T, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df,
         clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()
pdf(file="Final/Liver/CountMatrix_Liver_Inter.pdf", width=6, height=9)
pheatmap(assay(rld_liver.Sig3), cluster_rows=T, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df,
         clustering_distance_rows="euclidean", clustering_distance_cols="euclidean")
dev.off()




plot.all.data <- plotPCA(rld_liver.Sig1, intgroup = c("Treatment", "Temp"), returnData=TRUE, ntop=num.genes1)
percentVar <- round(100*attr(plot.all.data, "percentVar"))
pdf(file="Final/Liver/PCA_DGE_Temp_Liver.pdf", width=9, height=6)
ggplot(plot.all.data, aes(PC1, PC2, color=Treatment,shape=Temp)) +
  geom_point(size=4) +
  geom_text_repel(aes(label =name),
                  size = 4, max.overlaps = 32)+
  scale_color_manual(values=c(my.colors)) +
  theme_linedraw(base_size = 18) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
dev.off()

plot.all.data <- plotPCA(rld_liver.Sig2, intgroup = c("Treatment", "Temp"), returnData=TRUE, ntop=num.genes2)
percentVar <- round(100*attr(plot.all.data, "percentVar"))
pdf(file="Final/Liver/PCA_DGE_Treat_Liver.pdf", width=9, height=6)
ggplot(plot.all.data, aes(PC1, PC2, color=Treatment,shape=Temp)) +
  geom_point(size=4) +
  geom_text_repel(aes(label =name),
                  size = 4, max.overlaps = 32)+
  scale_color_manual(values=c(my.colors)) +
  theme_linedraw(base_size = 18) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
dev.off()

plot.all.data <- plotPCA(rld_liver.Sig3, intgroup = c("Treatment", "Temp"), returnData=TRUE, ntop=num.genes3)
percentVar <- round(100*attr(plot.all.data, "percentVar"))
pdf(file="Final/Liver/PCA_DGE_Inter_Liver.pdf", width=9, height=6)
ggplot(plot.all.data, aes(PC1, PC2, color=Treatment,shape=Temp)) +
  geom_point(size=4) +
  geom_text_repel(aes(label =name),
                  size = 4, max.overlaps = 32)+
  scale_color_manual(values=c(my.colors)) +
  theme_linedraw(base_size = 18) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
dev.off()

save.image("Final/Liver/through_logxform_liver.RData")
