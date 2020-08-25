library(ggplot2)
library(ggrepel)

#Set working directory
setwd("path/to/working/directory")

#get data
dat<-read.table("bactOG_counts_normalised.csv ", header=TRUE, sep=",", row.names=1, check.names=FALSE)

pca_results<- prcomp(dat, scale=TRUE)
#get just the PC numbers
pca_data<-as.data.frame(pca_results$rotation)
#force ordering so x-axis labels don't get moved around
pca_data$species<-factor(row.names(pca_data), levels=unique(row.names(pca_data)))
#add phenotype info to allow colouring later
pca_data$Phenotype<-c("HAP","HAP", "HAP", "SAP", "SAP", "SAP", "NAP", "NAP", "NAP")

#get info about the components: 
pc_summary<-summary(pca_results)

# Importance of components:
#   PC1    PC2     PC3     PC4     PC5     PC6    PC7     PC8     PC9
# Standard deviation     2.0551 1.0669 0.91754 0.79348 0.76222 0.71620 0.6722 0.58273 0.53062
# Proportion of Variance 0.4693 0.1265 0.09354 0.06996 0.06455 0.05699 0.0502 0.03773 0.03128
# Cumulative Proportion  0.4693 0.5957 0.68928 0.75923 0.82379 0.88078 0.9310 0.96872 1.00000

#get pc variance percentages:
pc1_var<-round((pc_summary$importance[2,1]) * 100, digits=1)
pc2_var<-round((pc_summary$importance[2,2]) * 100, digits=1)

#plot the pcs
pca_plot <- ggplot(pca_data, aes(PC1, PC2)) +
  geom_point(size=2, aes(color= Phenotype)) + 
  xlab(paste0("PC1: ",pc1_var,"% variance")) +
  ylab(paste0("PC2: ",pc2_var,"% variance")) +
  geom_text_repel(aes(label=species, color=Phenotype), size = 5, fontface="italic", point.padding=0.5, show.legend=FALSE) +
  theme(text=element_text(size=12 )) 


ggsave("HAP_NAP_SAP_COG_normalised_PCA.png", plot=pca_plot, device = "png", dpi = 300)


#look at the loadings
loadings = as.data.frame(pca_results$x)
compload = abs(loadings)
sweep(compload, 2, colSums(compload), "/")
#View(compload)
write.table(compload, file="component_loadings_scaledtosmallestgenome.table", sep="\t")