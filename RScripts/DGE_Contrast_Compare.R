###########################################################################
###                           Samarth Mathur, PhD                       ###
###                         The Ohio State University                   ###
###                                                                     ###
###     Date Created: 08/09/21                  Last Modified: 08/09/21 ###
###########################################################################
###########################################################################
###                   DGE_Interaction_Compare.R                         ###
###########################################################################

# Analysis based on : https://github.com/tavareshugo/tutorial_DESeq2_contrasts

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
library("BioVenn")

setwd("~/Documents/Postdoc_research/rnaseq/results/")

# Load data
genes <- read.csv("CroVir_GeneAnnotation.csv", header = T)
assembly <- readGFF("stringtie_all_merged.annotated.gtf")

inter_liver <- read.csv("Final/Liver/DGE/DGE_Con5_Liver.csv") # Temp x Treat interaction (Con vs Treat (at 26C in kidney)) N = 236
inter_kidney <- read.csv("Final/Kidney/DGE/DGE_Inter_Kidney.csv") # Temp x Treat interaction (Con vs Treat (at 26C in kidney)) N = 166


# Remove NAs
inter_liver <- inter_liver[-which(is.na(inter_liver$GeneID)),] # N = 91
inter_kidney <- inter_kidney[-which(is.na(inter_kidney$GeneID)),] # N = 67

# Get annotations for DEGs

inter_liver <- inner_join(genes, inter_liver, by=c("GeneID")) # 91
inter_kidney <- inner_join(genes, inter_kidney, by=c("GeneID")) # 67

inter_liver <- inter_liver[,-c(8,9)]
inter_kidney <- inter_kidney[,-c(8,9)]

write.csv(inter_liver,file="Final/Liver/genelist/Liver_DGE_Inter_withAnn.csv", quote = F, row.names = F)
write.csv(inter_kidney,file="Final/Kidney/genelist/Kidney_DGE_Inter_withAnn.csv", quote = F, row.names = F)

# Map changes due to interaction
liver_change <- c(dim(inter_liver[which(inter_liver$log2FoldChange > 0),])[1], # 41
                  dim(inter_liver[which(inter_liver$log2FoldChange < 0),])[1]) # 50

kidney_change <- c(dim(inter_kidney[which(inter_kidney$log2FoldChange > 0),])[1], # 56
                  dim(inter_kidney[which(inter_kidney$log2FoldChange < 0),])[1]) # 11

Regulation <- c("Up","Down")
cols <- c("#98BF64","#F69ABF")

inter_dge <- as.data.frame(cbind(c("Liver","Liver","Kidney","Kidney"),
                                 c("Up","Down","Up","Down"),
                                 c(liver_change,kidney_change)))

colnames(inter_dge) <- c("Tissue","Regulation","DEGs")
inter_dge$Regulation <- factor(inter_dge$Regulation,levels = c("Up", "Down"))
inter_dge$Tissue <- factor(inter_dge$Tissue,levels = c("Liver", "Kidney"))
inter_dge$DEGs <- as.numeric(inter_dge$DEGs)

pdf(file = "Tables and Figures/Interaction_DEGs.pdf", width = 7, height = 6)
ggplot(inter_dge, aes(fill=Regulation, y=DEGs, x=Tissue)) + 
  geom_bar(position="dodge", stat="identity",colour="black")+ylim(0,75)+
  theme_classic(base_size = 15)+scale_fill_manual(values=cols)+
  labs(fill = "Gene Regulation", x="Tissue", y="No. of differentially expressed genes (DEGs)",
       title = "Temperature X Treatment")
dev.off()


