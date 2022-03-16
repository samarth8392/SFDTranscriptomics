###########################################################################
###                           Samarth Mathur, PhD                       ###
###                         The Ohio State University                   ###
###                                                                     ###
###     Date Created: 07/14/21                  Last Modified: 07/30/21 ###
###########################################################################
###########################################################################
###                   gene_list.R   		                                ###
###########################################################################

setwd("~/Documents/Postdoc_research/rnaseq/results/Final/")
genes <- read.csv("../CroVir_GeneAnnotation.csv",header = T)

# Get genelist for each DEG (due to treatment only)
liver_dge_ann <- read.csv("Liver/genelist/Liver_DGE_treatonly_withAnn.csv") # 410
kidney_dge_ann <- read.csv("Kidney/genelist/Kidney_DGE_treatonly_withAnn.csv") # 415 
skin_dge_ann <- read.csv("Skin/genelist/Skin_DGE_treatonly_withAnn.csv") # 15

liver_dge_list <- unique(liver_dge_ann$Homolog_Anolis) # 337
kidney_dge_list <- unique(kidney_dge_ann$Homolog_Anolis) # 347
skin_dge_list <- unique(skin_dge_ann$Homolog_Anolis) #13

liver_dge_up <- liver_dge_ann[which(liver_dge_ann$log2FoldChange > 0),] # 268 
liver_dge_down <- liver_dge_ann[which(liver_dge_ann$log2FoldChange < 0),] # 142

kidney_dge_up <- kidney_dge_ann[which(kidney_dge_ann$log2FoldChange > 0),] # 130 
kidney_dge_down <- kidney_dge_ann[which(kidney_dge_ann$log2FoldChange < 0),] # 285

skin_dge_up <- skin_dge_ann[which(skin_dge_ann$log2FoldChange > 0),] # 12 
skin_dge_down <- skin_dge_ann[which(skin_dge_ann$log2FoldChange < 0),] # 3


write.csv(liver_dge_list, "Liver/genelist/Liver_DGE_treatonly_genelist.csv", quote = F, row.names = F, col.names = F)
write.csv(kidney_dge_list, "Kidney/genelist/Kidney_DGE_treatonly_genelist.csv", quote = F, row.names = F, col.names = F)
write.csv(skin_dge_list, "Skin/genelist/Skin_DGE_treatonly_genelist.csv", quote = F, row.names = F, col.names = F)

write.csv(unique(liver_dge_up$Homolog_Anolis), "Liver/genelist/Liver_DGE_treatonly_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(liver_dge_down$Homolog_Anolis), "Liver/genelist/Liver_DGE_treatonly_down_genelist.csv", quote = F, row.names = F)

write.csv(unique(kidney_dge_up$Homolog_Anolis), "Kidney/genelist/Kidney_DGE_treatonly_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(kidney_dge_down$Homolog_Anolis), "Kidney/genelist/Kidney_DGE_treatonly_down_genelist.csv", quote = F, row.names = F)

write.csv(unique(skin_dge_up$Homolog_Anolis), "Skin/genelist/Skin_DGE_treatonly_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(skin_dge_down$Homolog_Anolis), "Skin/genelist/Skin_DGE_treatonly_down_genelist.csv", quote = F, row.names = F)


# Get genelist for each DEG (due to temp only)

liver_temp_ann <- read.csv("Liver/genelist/Liver_DGE_TempOnly_withAnn.csv") # 134
kidney_temp_ann <- read.csv("Kidney/genelist/Kidney_DGE_TempOnly_withAnn.csv") # 95 

liver_temp_list <- unique(liver_temp_ann$Homolog_Anolis) # 118
kidney_temp_list <- unique(kidney_temp_ann$Homolog_Anolis) # 87

liver_temp_up <- liver_temp_ann[which(liver_temp_ann$log2FoldChange > 0),] # 70 
liver_temp_down <- liver_temp_ann[which(liver_temp_ann$log2FoldChange < 0),] # 64

kidney_temp_up <- kidney_temp_ann[which(kidney_temp_ann$log2FoldChange > 0),] # 38 
kidney_temp_down <- kidney_temp_ann[which(kidney_temp_ann$log2FoldChange < 0),] # 57

write.csv(unique(liver_temp_up$Homolog_Anolis), "Liver/genelist/Liver_DGE_TempOnly_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(liver_temp_down$Homolog_Anolis), "Liver/genelist/Liver_DGE_TempOnly_down_genelist.csv", quote = F, row.names = F)

write.csv(unique(kidney_temp_up$Homolog_Anolis), "Kidney/genelist/Kidney_DGE_TempOnly_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(kidney_temp_down$Homolog_Anolis), "Kidney/genelist/Kidney_DGE_TempOnly_down_genelist.csv", quote = F, row.names = F)

# Get genelist for each DEG (due to Interaction)

liver_inter_ann <- read.csv("Liver/genelist/Liver_DGE_Inter_withAnn.csv") # 91
kidney_inter_ann <- read.csv("Kidney/genelist/Kidney_DGE_Inter_withAnn.csv") # 67 

liver_inter_list <- unique(liver_inter_ann$Homolog_Anolis) # 73
kidney_inter_list <- unique(kidney_inter_ann$Homolog_Anolis) # 51

liver_inter_up <- liver_inter_ann[which(liver_inter_ann$log2FoldChange > 0),] # 50 
liver_inter_down <- liver_inter_ann[which(liver_inter_ann$log2FoldChange < 0),] # 41

kidney_inter_up <- kidney_inter_ann[which(kidney_inter_ann$log2FoldChange > 0),] # 56 
kidney_inter_down <- kidney_inter_ann[which(kidney_inter_ann$log2FoldChange < 0),] # 11

write.csv(unique(liver_inter_up$Homolog_Anolis), "Liver/genelist/Liver_DGE_Inter_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(liver_inter_down$Homolog_Anolis), "Liver/genelist/Liver_DGE_Inter_down_genelist.csv", quote = F, row.names = F)

write.csv(unique(kidney_inter_up$Homolog_Anolis), "Kidney/genelist/Kidney_DGE_Inter_up_genelist.csv", quote = F, row.names = F)
write.csv(unique(kidney_inter_down$Homolog_Anolis), "Kidney/genelist/Kidney_DGE_Inter_down_genelist.csv", quote = F, row.names = F)



# Get genelist for all modules

liver_mod_ann <- read.csv("Liver/genelist/Liver_Sigmods_withAnn.csv") # 1157
kidney_mod_ann <- read.csv("Kidney/genelist/Kidney_Sigmods_withAnn.csv") # 1532
skin_mod_ann <- read.csv("Skin/genelist/Skin_Sigmods_withAnn.csv") # 2061

liver_mod_list <- unique(liver_mod_ann$Homolog_Anolis) # 939
kidney_mod_list <- unique(kidney_mod_ann$Homolog_Anolis) # 1278
skin_mod_list <- unique(skin_mod_ann$Homolog_Anolis) #1658

write.csv(liver_mod_list, "Liver/genelist/Liver_Sigmods_genelist.csv", quote = F, row.names = F)
write.csv(kidney_mod_list, "Kidney/genelist/Kidney_Sigmods_genelist.csv", quote = F, row.names = F)
write.csv(skin_mod_list, "Skin/genelist/Skin_Sigmods_genelist.csv", quote = F, row.names = F)

# Get genelist for each modules

livermods <- unique(liver_mod_ann$Module)
kidneymods <- unique(kidney_mod_ann$Module)
skinmods <- unique(skin_mod_ann$Module)

for (i in skinmods)
{
  skin_mod_list <- unique(skin_mod_ann$Homolog_Anolis[which(skin_mod_ann$Module == i)])
  skin_mod_list <- skin_mod_list[-which(is.na(skin_mod_list))]
  write.csv(skin_mod_list, paste("Skin/genelist/Skin_",i,"_genelist.csv"), quote = F, row.names = F)
}


# Get genelist for each DGE + WGCNA

liver_both_ann <- read.csv("Liver/Liver_DGE_WGCNA_withAnn.csv") # 85
kidney_both_ann <- read.csv("Kidney/Kidney_DGE_WGCNA_withAnn.csv") # 105
skin_both_ann <- read.csv("Skin/Skin_DGE_WGCNA_withAnn.csv") # 1

liver_both_list <- unique(liver_both_ann$Homolog_Anolis) # 46
kidney_both_list <- unique(kidney_both_ann$Homolog_Anolis) # 68
skin_both_list <- unique(skin_both_ann$Homolog_Anolis) # 1

write.csv(liver_both_list, "Liver/Liver_both_genelist.csv", quote = F, row.names = F)
write.csv(kidney_both_list, "Kidney/Kidney_both_genelist.csv", quote = F, row.names = F)
write.csv(skin_both_list, "Skin/Skin_both_genelist.csv", quote = F, row.names = F)

