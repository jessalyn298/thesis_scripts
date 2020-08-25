library(ggplot2)
library(ggrepel)

setwd("path/to/working/directory")

# ###PCA1 on raw count data:###
# dat_raw<-read.delim("HAP_NAP_SAP_EGN_rawcounts.table.csv", sep=",", row.names=1)
# 
# pca_results<-prcomp(dat_raw, scale=TRUE)
# #get just the PC numbers
# pca_data<-as.data.frame(pca_results$rotation)
# #force ordering so x-axis labels don't get moved around
# pca_data$species<-factor(row.names(pca_data), levels=unique(row.names(pca_data)))
# #add phenotype info to allow colouring later
# pca_data$phenotype<-c("HAP","HAP", "HAP", "SAP", "SAP", "SAP", "NAP", "NAP", "NAP")
# 
# pca_summary<-(summary(pca_results)) #PRCOMP SCALING ON
# # Importance of components:
# #   PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8
# # Standard deviation     2.2892 1.0033 0.85357 0.77935 0.71208 0.54957 0.5161 0.47950
# # Proportion of Variance 0.5823 0.1119 0.08095 0.06749 0.05634 0.03356 0.0296 0.02555
# # Cumulative Proportion  0.5823 0.6941 0.77505 0.84254 0.89888 0.93244 0.9620 0.98758
# # PC9
# # Standard deviation     0.33431
# # Proportion of Variance 0.01242
# # Cumulative Proportion  1.00000
# 
# #get the variation for the PCs
# pc1_var<-round((pca_summary$importance[2,1]) * 100, digits=1)
# pc2_var<-round((pca_summary$importance[2,2]) *100, digits=1)
# #plot the PCA
# pca1_rawdata<-ggplot(pca_data, aes(PC1, PC2, color=phenotype)) +
#   geom_point(size=5) + 
#   xlab(paste0("PC1: ",pc1_var,"% variance")) +
#   ylab(paste0("PC2: ",pc2_var,"% variance")) +
#   # coord_fixed() +
#   geom_text_repel(aes(label=species)) +
#   ggtitle("PCA1 - raw counts")
# 
# #look at the loadings
# loadings = as.data.frame(pca_results$x)
# compload = abs(loadings)
# sweep(compload, 2, colSums(compload), "/")
# View(compload)
# #write.table(compload, file="component_loadings.table", sep="\t")
# 
# 
# 
# 
# ###PCA2 - scaled by genome size ###
# dat_gs<-read.delim("HAP_NAP_SAP_EGN_all_genomenormalised_counts.csv", sep=",", row.names=1)
# 
# pca_results_gs<-prcomp(dat_gs, scale=TRUE)
# #get just the PC numbers
# pca_data_gs<-as.data.frame(pca_results_gs$rotation)
# #force ordering so x-axis labels don't get moved around
# pca_data_gs$species<-factor(row.names(pca_data_gs), levels=unique(row.names(pca_data_gs)))
# #add phenotype info to allow colouring later
# pca_data_gs$phenotype<-c("HAP","HAP", "HAP", "SAP", "SAP", "SAP", "NAP", "NAP", "NAP")
# 
# pca_summary_gs<-(summary(pca_results_gs))
# # Importance of components:
# #   PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8
# # Standard deviation     2.2892 1.0033 0.85357 0.77935 0.71208 0.54957 0.5161 0.47950
# # Proportion of Variance 0.5823 0.1119 0.08095 0.06749 0.05634 0.03356 0.0296 0.02555
# # Cumulative Proportion  0.5823 0.6941 0.77505 0.84254 0.89888 0.93244 0.9620 0.98758
# # PC9
# # Standard deviation     0.33431
# # Proportion of Variance 0.01242
# # Cumulative Proportion  1.00000
# 
# #get the variation for the PCs
# pc1_var_gs<-round((pca_summary_gs$importance[2,1]) * 100, digits=1)
# pc2_var_gs<-round((pca_summary_gs$importance[2,2]) *100, digits=1)
# #plot the PCA
# pca2_gs<-ggplot(pca_data_gs, aes(PC1, PC2, color=phenotype)) +
#   geom_point(size=5) + 
#   xlab(paste0("PC1: ",pc1_var_gs,"% variance")) +
#   ylab(paste0("PC2: ",pc2_var_gs,"% variance")) +
#   # coord_fixed() +
#   geom_text_repel(aes(label=species))  
#   #ggtitle("PCA2 - scaled to genome size")
#   
# 
# #look at the loadings
# loadings_gs <- as.data.frame(pca_results_gs$x)
# compload_gs <- abs(loadings_gs)
# sweep(compload_gs, 2, colSums(compload_gs), "/")
# View(compload_gs)
# write.table(compload_gs, file="component_loadings.table", sep="\t")




###PCA3 - scaled using smallest genome###
dat_smolg<-read.delim("HAP_NAP_SAP_EGN_counts_scaledtosmallestgenome.csv", sep=",", row.names=1, check.names = FALSE)

pca_results_smolg<-prcomp(dat_smolg, scale=TRUE)
#get just the PC numbers
pca_data_smolg<-as.data.frame(pca_results_smolg$rotation)
#force ordering so x-axis labels don't get moved around
pca_data_smolg$species<-factor(row.names(pca_data_smolg), levels=unique(row.names(pca_data_smolg)))
#add phenotype info to allow colouring later
pca_data_smolg$Phenotype<-c("HAP","HAP", "HAP", "SAP", "SAP", "SAP", "NAP", "NAP", "NAP")

pca_summary_smolg<-(summary(pca_results_smolg))
# Importance of components:
#   PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9
# Standard deviation     2.2892 1.0033 0.85357 0.77935 0.71208 0.54957 0.5161 0.47950 0.33431
# Proportion of Variance 0.5823 0.1119 0.08095 0.06749 0.05634 0.03356 0.0296 0.02555 0.01242
# Cumulative Proportion  0.5823 0.6941 0.77505 0.84254 0.89888 0.93244 0.9620 0.98758 1.00000

#get the variation for the PCs
pc1_var_smolg<-round((pca_summary_smolg$importance[2,1]) * 100, digits=1)
pc2_var_smolg<-round((pca_summary_smolg$importance[2,2]) *100, digits=1)
#plot the PCA
pca3_smolg<- ggplot(pca_data_smolg, aes(PC1, PC2)) +
  geom_point(size=2, aes(color= Phenotype)) + 
  xlab(paste0("PC1: ",pc1_var_smolg,"% variance")) +
  ylab(paste0("PC2: ",pc2_var_smolg,"% variance")) +
  geom_text_repel(aes(label=species, color=Phenotype), size = 5, fontface="italic", point.padding=0.5, show.legend=FALSE) +
  theme(text=element_text(size=12 )) 
 

ggsave(filename = "EGN_allcounts_scaledtosmallestgenome_PCAplot.png", device="png", dpi = 300)


#look at the loadings
loadings_smolg = as.data.frame(pca_results_smolg$x)
compload_smolg = abs(loadings_smolg)
sweep(compload_smolg, 2, colSums(compload_smolg), "/")
View(compload_smolg)
write.table(compload_smolg, file="component_loadings_scaledtosmallestgenome.table", sep="\t")
