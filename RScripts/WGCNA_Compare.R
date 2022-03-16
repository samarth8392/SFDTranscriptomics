###########################################################################
###                           Samarth Mathur, PhD                       ###
###                         The Ohio State University                   ###
###                                                                     ###
###     Date Created: 06/29/21                  Last Modified: 07/14/21 ###
###########################################################################
###########################################################################
###                   WGCNA_Compare.R   		                            ###
###########################################################################

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
library("gridExtra")
library("WGCNA")
library(ggpubr)
setwd("~/Documents/Postdoc_research/rnaseq/results/Final/")

annotation <- read.csv("../Genenames.csv", header = T)
genes <- read.csv("../CroVir_GeneAnnotation.csv",header = T)

liver_mods <- read.csv("Liver/WGCNA/Liver_Module_gene_IDs_noNA.csv") # N =1157
kidney_mods <- read.csv("Kidney/WGCNA/Kidney_Module_gene_IDs_noNA.csv") # N = 1532
skin_mods <- read.csv("Skin/WGCNA/Skin_Module_gene_IDs_noNA.csv") # N = 2061

liver_dge <- read.csv("Liver/DGEs_TreatOnly_Liver.csv") # N = 410
kidney_dge <- read.csv("Kidney/DGEs_TreatOnly_Kidney.csv") # N =415
skin_dge <- read.csv("Skin/DGEs_TreatOnly_Skin.csv") # N =15

liver_both <- inner_join(liver_mods,liver_dge, by="GeneID") # N =85
kidney_both <- inner_join(kidney_mods,kidney_dge, by="GeneID") # N = 105
skin_both <- inner_join(skin_mods,skin_dge, by="GeneID") # N = 1

liver_both <- liver_both[,-c(3,4,7,8)]
kidney_both <- kidney_both[,-c(3,4,7,8)]
skin_both <- skin_both[,-c(3,4,7,8)]

liver_both$Module <- factor(liver_both$Module, levels=c("darkgreen","burlywood","antiquewhite1","magenta2"))
kidney_both$Module <- factor(kidney_both$Module, levels=c("magenta","navajowhite","magenta4","green2"))

# Plot DGEs within different modules

p1a <- ggplot(liver_both)+ylim(-5,5)+
  geom_boxplot(aes(x=Module, y=log2FoldChange, fill=Module), size=1)+
  theme_classic(base_size = 15)+
  labs(x="Module", y =expression('log'['2']~"(fold change)"))+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray50", size=1.5)+
  scale_fill_manual(values=c("darkgreen","burlywood", "antiquewhite1","magenta2"))

p1b <- ggplot(kidney_both)+ylim(-5,5)+
  geom_boxplot(aes(x=Module, y=log2FoldChange, fill=Module), size=1)+
  theme_classic(base_size = 15)+
  labs(x="Module", y =expression('log'['2']~"(fold change)"))+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray50", size=1.5)+
  scale_fill_manual(values=c("magenta","navajowhite","magenta4","green2"))

p1c <- ggplot(skin_both)+ylim(-5,5)+
  geom_boxplot(aes(x=Module, y=log2FoldChange, fill=Module), size=1)+
  theme_classic(base_size = 15)+
  labs(x="Module", y =expression('log'['2']~"(fold change)"))+
  theme(legend.position = "none")+
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray50", size=1.5)

grid.arrange(p1a, p1b, p1c, nrow = 1, widths = c(1.5,1.5,0.6))


p2a <- ggplot(liver_both, aes(x=Gene_MM, y=Gene_TS, color=Module)) +
  geom_point(size=3, shape=16) +
  geom_point(shape = 1,size=3,colour = "black") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(x="Module Membership", y ="Trait Significance",title = "Liver")+
  theme(legend.position = "none")+
  theme_classic(base_size = 15)+
  scale_color_manual(values=c("darkgreen","burlywood", "antiquewhite1","magenta2"))
p2a

summary(lm(liver_both$Gene_TS[which(liver_both$Module=="darkgreen")] ~ 
             liver_both$Gene_MM[which(liver_both$Module=="darkgreen")])) # p-value: 0.00701, R2 = 0.2891

summary(lm(liver_both$Gene_TS[which(liver_both$Module=="burlywood")] ~ 
             liver_both$Gene_MM[which(liver_both$Module=="burlywood")])) # p-value: 0.1009, R2 = 0.0587 

summary(lm(liver_both$Gene_TS[which(liver_both$Module=="antiquewhite1")] ~ 
             liver_both$Gene_MM[which(liver_both$Module=="antiquewhite1")]))  # p-value: 9.753e-08, R2 = 0.7202 

summary(lm(liver_both$Gene_TS[which(liver_both$Module=="magenta2")] ~ 
             liver_both$Gene_MM[which(liver_both$Module=="magenta2")])) # p-value: 0.1873, R2 = 0.1243 

p2b <- ggplot(kidney_both, aes(x=Gene_MM, y=Gene_TS, color=Module)) +
  geom_point(size=3, shape=16) +
  geom_point(shape = 1,size=3,colour = "black") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(x="Module Membership", y ="Trait Significance",title = "Kidney")+
  theme(legend.position = "none")+
  theme_classic(base_size = 15)+
  scale_color_manual(values=c("magenta","navajowhite","magenta4","green2"))
p2b

summary(lm(kidney_both$Gene_TS[which(kidney_both$Module=="magenta")] ~ 
             kidney_both$Gene_MM[which(kidney_both$Module=="magenta")])) # p-value: 0.0624, R2 = 0.02897

summary(lm(kidney_both$Gene_TS[which(kidney_both$Module=="navajowhite")] ~ 
             kidney_both$Gene_MM[which(kidney_both$Module=="navajowhite")])) # p-value: 0.1009, R2 = 0.0587 

summary(lm(kidney_both$Gene_TS[which(kidney_both$Module=="magenta4")] ~ 
             kidney_both$Gene_MM[which(kidney_both$Module=="magenta4")]))  # p-value: NA, R2 = NA 

summary(lm(kidney_both$Gene_TS[which(kidney_both$Module=="green2")] ~ 
             kidney_both$Gene_MM[which(kidney_both$Module=="green2")])) # p-value: 0.2948, R2 = 0.02773 


p2c <- ggplot(skin_both, aes(x=Gene_MM, y=Gene_TS, color=Module)) +
  geom_point(size=3, shape=16) +
  geom_point(shape = 1,size=3,colour = "black") +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+
  labs(x="Module Membership", y ="Trait Significance",title = "Kidney")+
  theme(legend.position = "none")+
  theme_classic(base_size = 15)+
  scale_color_manual(values=c("firebrick4"))
p2c


# Get annotation for each DEG
liver_dge_ann <- inner_join(genes, liver_dge, by="GeneID")
kidney_dge_ann <- inner_join(genes, kidney_dge, by="GeneID")
skin_dge_ann <- inner_join(genes, skin_dge, by="GeneID")

liver_dge_ann <- liver_dge_ann[,-c(8,9)]
kidney_dge_ann <- kidney_dge_ann[,-c(8,9)]
skin_dge_ann <- skin_dge_ann[,-c(8,9)]
write.csv(liver_dge_ann, "Liver/Liver_DGE_treatonly_withAnn.csv", quote = F, row.names = F)
write.csv(kidney_dge_ann, "Kidney/Kidney_DGE_treatonly_withAnn.csv", quote = F, row.names = F)
write.csv(skin_dge_ann, "Skin/Skin_DGE_treatonly_withAnn.csv", quote = F, row.names = F)

# Get annotation for each module
liver_mod_ann <- inner_join(genes, liver_mods, by="GeneID")
kidney_mod_ann <- inner_join(genes, kidney_mods, by="GeneID")
skin_mod_ann <- inner_join(genes, skin_mods, by="GeneID")

liver_mod_ann <- liver_mod_ann[,-c(9,10)]
kidney_mod_ann <- kidney_mod_ann[,-c(9,10)]
skin_mod_ann <- skin_mod_ann[,-c(9,10)]
write.csv(liver_mod_ann, "Liver/Liver_Sigmods_withAnn.csv", quote = F, row.names = F)
write.csv(kidney_mod_ann, "Kidney/Kidney_Sigmods_withAnn.csv", quote = F, row.names = F)
write.csv(skin_mod_ann, "Skin/Skin_Sigmods_withAnn.csv", quote = F, row.names = F)

# Get annotation for DGE and module
liver_both_ann <- inner_join(genes, liver_both, by="GeneID")
kidney_both_ann <- inner_join(genes, kidney_both, by="GeneID")
skin_both_ann <- inner_join(genes, skin_both, by="GeneID")

write.csv(liver_both_ann, "Liver/Liver_DGE_WGCNA_withAnn.csv", quote = F, row.names = F)
write.csv(kidney_both_ann, "Kidney/Kidney_DGE_WGCNA_withAnn.csv", quote = F, row.names = F)
write.csv(skin_both_ann, "Skin/Skin_DGE_WGCNA_withAnn.csv", quote = F, row.names = F)



