###########################################################################
###                           Samarth Mathur, PhD                       ###
###                         The Ohio State University                   ###
###                                                                     ###
###     Date Created: 08/09/21                  Last Modified: 08/09/21 ###
###########################################################################
###########################################################################
###                   DGE_Compare_Temp.R                                ###
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
library("BioVenn")

setwd("~/Documents/Postdoc_research/rnaseq/results/")

# Load data
temp_liver <- read.csv("Final/Liver/DGE/DGE_Temp_Liver.csv")
treat_liver <- read.csv("Final/Liver/DGE/DGE_Treat_Liver.csv")
inter_liver <- read.csv("Final/Liver/DGE/DGE_Inter_Liver.csv")

temp_kidney <- read.csv("Final/Kidney/DGE/DGE_Temp_Kidney.csv")
treat_kidney <- read.csv("Final/Kidney/DGE/DGE_Treat_Kidney.csv")
inter_kidney <- read.csv("Final/Kidney/DGE/DGE_Inter_Kidney.csv")

temp_skin <- read.csv("Final/skin/DGE/DGE_Temp_skin.csv")
treat_skin <- read.csv("Final/skin/DGE/DGE_Treat_skin.csv")

# Remove NAs
temp_liver <- temp_liver[-which(is.na(temp_liver$GeneID)),] # N = 182
treat_liver <- treat_liver[-which(is.na(treat_liver$GeneID)),] # N =503
inter_liver <- inter_liver[-which(is.na(inter_liver$GeneID)),] # N = 91

temp_kidney <- temp_kidney[-which(is.na(temp_kidney$GeneID)),] # N = 131
treat_kidney <- treat_kidney[-which(is.na(treat_kidney$GeneID)),] # N =507
inter_kidney <- inter_kidney[-which(is.na(inter_kidney$GeneID)),] # N = 67

treat_skin <- treat_skin[-which(is.na(treat_skin$GeneID)),] # N=15

# Partition data
# Liver
df1 <- inner_join(temp_liver, treat_liver, by=c("GeneID"))
df2 <- inner_join(treat_liver, inter_liver, by=c("GeneID"))
df3 <- inner_join(temp_liver, inter_liver, by=c("GeneID"))
df4 <- inner_join(df1, df2, df3, by=c("GeneID"))
df1 <- anti_join(df1,df4,by=c("GeneID")) # N = 19
df2 <- anti_join(df2,df4,by=c("GeneID")) # N = 55
df3 <- anti_join(df3,df4,by=c("GeneID")) # N = 7

df5 <- anti_join(treat_liver,df1,by=c("GeneID")) # Treat - Temp+treat (N=484)
df6 <- anti_join(df5,df2,by=c("GeneID")) # Treat - Temp+Treat - Treat+Inter (N=433)
treat_only_liver <- anti_join(df6,df4,by=c("GeneID")) #  Treat - Temp+Treat - Treat+Inter - All (N=410)

df7 <- anti_join(temp_liver,df1,by=c("GeneID")) # Temp - Temp+Treat (N=164)
df8 <- anti_join(df7,df3,by=c("GeneID")) # Temp - Temp+Treat - Temp+Inter (N=159)
temp_only_liver <- anti_join(df8,df4,by=c("GeneID")) # Temp - Temp+Treat - Temp+Inter - All (N=134)

df9 <- anti_join(inter_liver,df2,by=c("GeneID")) # Inter - Treat+Inter (N=43)
df10 <- anti_join(df9,df3,by=c("GeneID")) # Inter - Treat+Inter - Temp+Inter (N=38)
inter_only_liver <- anti_join(df10,df4,by=c("GeneID")) # Inter - Treat+Inter - Temp+Inter - All (N=14)


# Kidney
df1 <- inner_join(temp_kidney, treat_kidney, by=c("GeneID"))
df2 <- inner_join(treat_kidney, inter_kidney, by=c("GeneID"))
df3 <- inner_join(temp_kidney, inter_kidney, by=c("GeneID"))
df4 <- inner_join(df1, df2, df3, by=c("GeneID"))
df1 <- anti_join(df1,df4,by=c("GeneID")) # N = 27
df2 <- anti_join(df2,df4,by=c("GeneID")) # N = 58
df3 <- anti_join(df3,df4,by=c("GeneID")) # N = 2

df5 <- anti_join(treat_kidney,df1,by=c("GeneID")) # Treat - Temp+treat (N=481)
df6 <- anti_join(df5,df2,by=c("GeneID")) # Treat - Temp+Treat - Treat+Inter (N=425)
treat_only_kidney <- anti_join(df6,df4,by=c("GeneID")) #  Treat - Temp+Treat - Treat+Inter - All (N=415)

df7 <- anti_join(temp_kidney,df1,by=c("GeneID")) # Temp - Temp+Treat (N=105)
df8 <- anti_join(df7,df3,by=c("GeneID")) # Temp - Temp+Treat - Temp+Inter (N=103)
temp_only_kidney <- anti_join(df8,df4,by=c("GeneID")) # Temp - Temp+Treat - Temp+Inter - All (N=95)

df9 <- anti_join(inter_kidney,df2,by=c("GeneID")) # Inter - Treat+Inter (N=19)
df10 <- anti_join(df9,df3,by=c("GeneID")) # Inter - Treat+Inter - Temp+Inter (N=17)
inter_only_kidney <- anti_join(df10,df4,by=c("GeneID")) # Inter - Treat+Inter - Temp+Inter - All (N=8)


#skin
treat_only_skin <- anti_join(treat_skin, temp_skin) # N =15

draw.venn(temp_liver$StringID,treat_liver$StringID,inter_liver$StringID,
          title="DGEs in Liver", subtitle = NULL,
          xtitle="Temp", ytitle="Treat", ztitle = "Inter")

draw.venn(temp_kidney$StringID,treat_kidney$StringID,inter_kidney$StringID,
          title="DGEs in Kidney", subtitle = NULL,
          xtitle="Temp", ytitle="Treat", ztitle = "Inter")

draw.venn(temp_skin$StringID,treat_skin$StringID,NULL,
          title="DGEs in Kidney", subtitle = NULL,
          xtitle="Temp", ytitle="Treat", ztitle = NULL)

draw.venn(treat_only_liver$StringID,treat_only_kidney$StringID,treat_only_skin$StringID,
          title="DGEs (Treatment Only)", subtitle = NULL,
          xtitle="Liver", ytitle="Kidney", ztitle = "Skin")


# Map changes due to fixed effects

condition <- c("Temperature","Treatment", "Temperature x Treatment")
change <- c("Up","Down")
cols <- c("#98BF64","#F69ABF")
liver_change <- c(dim(temp_liver[which(temp_liver$log2FoldChange > 0),])[1],
                  dim(treat_liver[which(treat_liver$log2FoldChange < 0),])[1],
                  dim(inter_liver[which(inter_liver$log2FoldChange > 0),])[1],
                  dim(temp_liver[which(temp_liver$log2FoldChange < 0),])[1],
                  dim(treat_liver[which(treat_liver$log2FoldChange > 0),])[1],
                  dim(inter_liver[which(inter_liver$log2FoldChange < 0),])[1])
liver_dge <- data.frame(condition,change,liver_change)
liver_dge$condition <- factor(liver_dge$condition,levels = c("Treatment", "Temperature", "Temperature x Treatment"))
liver_dge$change <- factor(liver_dge$change,levels = c("Up", "Down"))

kidney_change <- c(dim(temp_kidney[which(temp_kidney$log2FoldChange > 0),])[1],
                   dim(treat_kidney[which(treat_kidney$log2FoldChange < 0),])[1],
                   dim(inter_kidney[which(inter_kidney$log2FoldChange > 0),])[1],
                   dim(temp_kidney[which(temp_kidney$log2FoldChange < 0),])[1],
                   dim(treat_kidney[which(treat_kidney$log2FoldChange > 0),])[1],
                   dim(inter_kidney[which(inter_kidney$log2FoldChange < 0),])[1])
kidney_dge <- data.frame(condition,change,kidney_change)
kidney_dge$condition <- factor(kidney_dge$condition,levels = c("Treatment", "Temperature", "Temperature x Treatment"))
kidney_dge$change <- factor(kidney_dge$change,levels = c("Up", "Down"))

condition2 <- c("Temperature","Treatment","Treatment","Temperature")
skin_change <- c(dim(temp_skin[which(temp_skin$log2FoldChange > 0),])[1],
                 dim(treat_skin[which(treat_skin$log2FoldChange < 0),])[1],
                 dim(treat_skin[which(treat_skin$log2FoldChange > 0),])[1],
                 dim(temp_skin[which(temp_skin$log2FoldChange < 0),])[1])
skin_dge <- data.frame(condition2,change,skin_change)
skin_dge$condition2 <- factor(skin_dge$condition2,levels = c("Treatment", "Temperature"))
skin_dge$change <- factor(skin_dge$change,levels = c("Up", "Down"))


pdf(file = "Final/Skin/DGEsSkin.pdf", width = 5, height = 6)
ggplot(skin_dge, aes(fill=change, y=skin_change, x=condition2)) + 
  geom_bar(position="dodge", stat="identity",colour="black")+ylim(0,15)+
  theme_classic(base_size = 15)+scale_fill_manual(values=cols)+
  labs(fill = "Gene Regulation", x="Fixed effect", y="No. of differentially expressed genes (DEGs)")
dev.off()

# Map changes due to fixed effects ONLY

condition <- c("Temperature","Treatment", "Temperature x Treatment")
change <- c("Up","Down")
cols <- c("#98BF64","#F69ABF")
liver_change <- c(dim(temp_only_liver[which(temp_only_liver$log2FoldChange > 0),])[1], # 70
                  dim(treat_only_liver[which(treat_only_liver$log2FoldChange < 0),])[1], # 142
                  dim(inter_only_liver[which(inter_only_liver$log2FoldChange > 0),])[1], # 11
                  dim(temp_only_liver[which(temp_only_liver$log2FoldChange < 0),])[1], # 64
                  dim(treat_only_liver[which(treat_only_liver$log2FoldChange > 0),])[1], # 268
                  dim(inter_only_liver[which(inter_only_liver$log2FoldChange < 0),])[1]) # 3
liver_dge <- data.frame(condition,change,liver_change)
liver_dge$condition <- factor(liver_dge$condition,levels = c("Treatment", "Temperature", "Temperature x Treatment"))
liver_dge$change <- factor(liver_dge$change,levels = c("Up", "Down"))

pdf(file = "Final/Liver/DGEs_Only_Liver.pdf", width = 8, height = 6)
ggplot(liver_dge, aes(fill=change, y=liver_change, x=condition)) + 
  geom_bar(position="dodge", stat="identity",colour="black")+ylim(0,300)+
  theme_classic(base_size = 15)+scale_fill_manual(values=cols)+
  labs(fill = "Gene Regulation", x="Fixed effect (Only)", y="No. of differentially expressed genes (DEGs)",
       title = "Liver")
dev.off()

kidney_change <- c(dim(temp_only_kidney[which(temp_only_kidney$log2FoldChange > 0),])[1], # 38
                  dim(treat_only_kidney[which(treat_only_kidney$log2FoldChange < 0),])[1], # 285
                  dim(inter_only_kidney[which(inter_only_kidney$log2FoldChange > 0),])[1], # 8
                  dim(temp_only_kidney[which(temp_only_kidney$log2FoldChange < 0),])[1], # 57
                  dim(treat_only_kidney[which(treat_only_kidney$log2FoldChange > 0),])[1], # 130
                  dim(inter_only_kidney[which(inter_only_kidney$log2FoldChange < 0),])[1]) # 0
kidney_dge <- data.frame(condition,change,kidney_change)
kidney_dge$condition <- factor(kidney_dge$condition,levels = c("Treatment", "Temperature", "Temperature x Treatment"))
kidney_dge$change <- factor(kidney_dge$change,levels = c("Up", "Down"))

pdf(file = "Final/Kidney/DGEs_Only_Kidney.pdf", width = 8, height = 6)
ggplot(kidney_dge, aes(fill=change, y=kidney_change, x=condition)) + 
  geom_bar(position="dodge", stat="identity",colour="black")+ylim(0,300)+
  theme_classic(base_size = 15)+scale_fill_manual(values=cols)+
  labs(fill = "Gene Regulation", x="Fixed effect (Only)", y="No. of differentially expressed genes (DEGs)",
       title = "Kidney")
dev.off()



condition2 <- c("Temperature","Treatment","Treatment","Temperature")
skin_change <- c(dim(temp_skin[which(temp_skin$log2FoldChange > 0),])[1], # 1
                 dim(treat_skin[which(treat_skin$log2FoldChange < 0),])[1], # 3
                 dim(treat_skin[which(treat_skin$log2FoldChange > 0),])[1], # 12
                 dim(temp_skin[which(temp_skin$log2FoldChange < 0),])[1]) # 1
skin_dge <- data.frame(condition2,change,skin_change)
skin_dge$condition2 <- factor(skin_dge$condition2,levels = c("Treatment", "Temperature"))
skin_dge$change <- factor(skin_dge$change,levels = c("Up", "Down"))

pdf(file = "Final/Skin/DGEsSkin.pdf", width = 5, height = 6)
ggplot(skin_dge, aes(fill=change, y=skin_change, x=condition2)) + 
  geom_bar(position="dodge", stat="identity",colour="black")+ylim(0,15)+
  theme_classic(base_size = 15)+scale_fill_manual(values=cols)+
  labs(fill = "Gene Regulation", x="Fixed effect (Only)", y="No. of differentially expressed genes (DEGs)",
       title = "Skin")
dev.off()

# Get annotations for DEGs

genes <- read.csv("CroVir_GeneAnnotation.csv",header = T)

#Liver
temp_only_liver <- inner_join(genes, temp_only_liver, by=c("GeneID")) # 134
treat_only_liver <- inner_join(genes, treat_only_liver, by=c("GeneID")) # 410
inter_only_liver <- inner_join(genes, inter_only_liver, by=c("GeneID")) # 14

temp_only_liver <- temp_only_liver[,-c(8,9)]
treat_only_liver <- treat_only_liver[,-c(8,9)]
inter_only_liver <- inter_only_liver[,-c(8,9)]

write.csv(temp_only_liver,file="Final/Liver/genelist/Liver_DGE_TempOnly_withAnn.csv", quote = F, row.names = F)
write.csv(treat_only_liver,file="Final/Liver/genelist/Liver_DGE_TreatOnly_withAnn.csv", quote = F, row.names = F)
write.csv(inter_only_liver,file="Final/Liver/genelist/Liver_DGE_InterOnly_withAnn.csv", quote = F, row.names = F)

#Kidney
temp_only_kidney <- inner_join(genes, temp_only_kidney, by=c("GeneID")) # 95
treat_only_kidney <- inner_join(genes, treat_only_kidney, by=c("GeneID")) # 415
inter_only_kidney <- inner_join(genes, inter_only_kidney, by=c("GeneID")) # 8

temp_only_kidney <- temp_only_kidney[,-c(8,9)]
treat_only_kidney <- treat_only_kidney[,-c(8,9)]
inter_only_kidney <- inter_only_kidney[,-c(8,9)]

write.csv(temp_only_kidney,file="Final/Kidney/genelist/Kidney_DGE_TempOnly_withAnn.csv", quote = F, row.names = F)
write.csv(treat_only_kidney,file="Final/Kidney/genelist/Kidney_DGE_TreatOnly_withAnn.csv", quote = F, row.names = F)
write.csv(inter_only_kidney,file="Final/Kidney/genelist/Kidney_DGE_InterOnly_withAnn.csv", quote = F, row.names = F)

#Skin
temp_skin <- inner_join(genes, temp_skin, by=c("GeneID")) # 0
treat_skin <- inner_join(genes, treat_skin, by=c("GeneID")) # 15

write.csv(treat_skin,file="Final/Skin/genelist/Skin_DGE_TreatOnly_withAnn.csv", quote = F, row.names = F)