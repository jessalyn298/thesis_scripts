library(ggplot2)
library(reshape2)
library(ggpubr)

# [1] ggpubr_0.2.5       reshape2_1.4.3          ggplot2_3.1.1         


setwd("path/to/working/directory")
dat <- read.delim("TPM_aatransporters.tsv", check.names=FALSE)


dat1<- melt(dat, id=c("sample", "phenotype"))
#This 'melts' the data, turning it into a table like this:
#bacteria 1 group 1 gene 1 3
#bacteria 1 group 1 gene 2 1
#bacteria 2 group 1 gene 1 1
#bacteria 2 group 1 gene 2 5
#etc

#to assign colours (I need 21 colours as I have 21 transporters)
mycols <-c("navy", "aquamarine", "deepskyblue4", "darkorchid4",
           "burlywood",   
           "firebrick1", "lightpink", "maroon4", "deeppink2", 
           "beige",
           "darkolivegreen3", "limegreen", "seagreen",
           "black") 



stack<-ggplot(dat1, aes(x = sample, y = value, fill = variable)) +
  geom_bar(stat = "identity")+
  facet_wrap(~phenotype, scales = "free_x") +
  scale_fill_manual(values = mycols)+
  guides(fill=guide_legend(title="Transporter Family")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        plot.margin = margin(10, 10, 10, 55),
        axis.title.y=element_text(margin=margin(t = 0, r = 10, b = 0, l = 0)),
        text = element_text(size=12)) +
  labs(x = "Species and Transcriptome Sample") +
  labs(y = "Transcripts per million total") 

ggsave(plot=stack, "TPM_aatransporters_stackedgraph_feb2020.png", device="png", dpi = 300)



####boxplots###

bp_dat<-read.delim("TPM_aatransporters_boxplots.tsv", check.names=FALSE)
bp_dat_melt<-melt(bp_dat, id=c("Species", "phenotype"))
x = as.factor(bp_dat$Species)

dummy <- read.delim("maximum_tpm_dummydata.txt", check.names=FALSE)
dummy$Species <- "A.sticklandii 12662"
dummy<-melt(dummy) # this makes dummy data which is the maximum point plotted per variable times 1.2 so 
                    # it can be plotted as an invisible point to increase the y axis. 


bp<-ggplot(bp_dat_melt, aes(x=variable, y=value, fill=phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(data=dummy, alpha=0.0, aes(fill=NA)) + #this adds an invisible point to increase the y-axis in all plots.
  stat_compare_means(aes(group = phenotype, label = paste0("P = ", ..p.format..)), method = "wilcox.test", 
                     hjust=0.5, size = 3.5) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#008000", "#3A9494","#80FF00")) +
  facet_wrap(~variable, scales='free', nrow = 4, ncol = 4) +
  theme(axis.text.x = element_blank(),
        text = element_text(size=12)) +
  labs(x = "Species") +
  labs(y = "Transcript per Million Total") 
  
ggsave(plot=bp, "TPM_boxplots_aatransporters.png", device="png", dpi = 300)
 
###Significance testing #### 
#first test normality:
list_of_transporters <- (colnames(bp_dat))[-c(1,2)]
for (i in list_of_transporters) {
  print(i)
  print(shapiro.test(as.matrix(bp_dat[i])))
}

# [1] "2.A.118."
# W = 0.46724, p-value = 2.752e-08
# [1] "2.A.3."
# W = 0.63576, p-value = 1.589e-06
# [1] "3.A.1.3."
# W = 0.82366, p-value = 0.0007368
# [1] "3.A.1.4."
# W = 0.74473, p-value = 4.176e-05
# [1] "2.A.120."
# W = 0.47779, p-value = 3.456e-08 
# [1] "2.A.23."
# W = 0.8107, p-value = 0.000443
# [1] "2.A.26."
# W = 0.47588, p-value = 3.315e-08 
# [1] "2.A.21.2."
# W = 0.62282, p-value = 1.122e-06
# [1] "2.A.25."
# W = 0.85064, p-value = 0.002247
# [1] "2.A.17."
# W = 0.48047, p-value = 3.663e-08
# [1] "2.A.78."
# W = 0.83877, p-value = 0.001362
# [1] "2.A.114."
# W = 0.82316, p-value = 0.0007221
# [1] "2.A.67."
# W = 0.47406, p-value = 3.187e-08
# [1] "1.A.11."
# W = 0.58507, p-value = 4.227e-07

#all are <0.05, implying the data is not normal.

#This way will give a letter to each species and test significance between them
# dat_12662 <- boxplot_dat[grep("12662", boxplot_dat$sample),]
# dat_49906 <- boxplot_dat[grep("49906", boxplot_dat$sample),]
# dat_isol6 <- boxplot_dat[grep("isol6", boxplot_dat$sample),]
# dat_S85 <- boxplot_dat[grep("S85", boxplot_dat$sample),]
# dat_SY3 <- boxplot_dat[grep("SY3", boxplot_dat$sample),]
# dat_007c <- boxplot_dat[grep("007c", boxplot_dat$sample),]
# 

dat_HAPs <- bp_dat[grep("HAP", bp_dat$phenotype),]
dat_NAPs <- bp_dat[grep("NAP", bp_dat$phenotype),]

list_of_transporters <- (colnames(boxplot_dat))[-c(1,2)]

for (i in list_of_transporters) {
  print(i)
  print(wilcox.test(as.matrix(dat_HAPs[i]),as.matrix(dat_NAPs[i]), paired = FALSE, alternative = "two.sided" ))
}

# [1] "2.A.118."
# W = 96, p-value = 0.03669
# [1] "2.A.3."
# W = 120, p-value = 0.001084
# [1] "3.A.1.3."
# W = 128, p-value = 0.000656
# [1] "3.A.1.4."
# W = 94, p-value = 0.1847
# [1] "2.A.120."
# W = 96, p-value = 0.03669
# [1] "2.A.23."
# W = 144, p-value = 1.027e-05
# [1] "2.A.26."
# W = 96, p-value = 0.03669
# [1] "2.A.21.2."
# W = 96, p-value = 0.1739
# [1] "2.A.25."
# W = 144, p-value = 7.396e-07
# [1] "2.A.17."
# W = 120, p-value = 0.001084
# [1] "2.A.78."
# W = 24, p-value = 0.005208
# [1] "2.A.114."
# W = 144, p-value = 2.604e-05
# [1] "2.A.67."
# W = 96, p-value = 0.03669
# [1] "1.A.11."
# W = 48, p-value = 0.1669


