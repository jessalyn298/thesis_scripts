setwd("C:/Users/Jess/OneDrive - Aberystwyth University/Lab/results/VFAs")
concs <- read.table("VFA_raw_concs.txt", sep="\t", header=TRUE, check.names=FALSE)

library(reshape2)
library(agricolae)
library(ggpubr)

# > sessionInfo()
# R version 3.6.0 (2019-04-26)
# [1] ggpubr_0.2.5    magrittr_1.5    ggplot2_3.1.1   agricolae_1.3-1 reshape2_1.4.3 

melt_concs <- melt(data = concs, id=c("Sample"))

split <- split(melt_concs, melt_concs$variable)
total<- split$Total
acetic <- split$`mM Acetic`
propanoic <-split$`mM Propionic`
butyric <- split$`mM Butyric total`
valeric <- split$`mM Valeric total`
caproic <- split$`mM Caproic total`

shapiro.test(total$value)
#W = 0.71261, p-value = 1.565e-09
shapiro.test(acetic$value)
#W = 0.91251, p-value = 0.000388
shapiro.test(propanoic$value)
#W = 0.66628, p-value = 2.095e-10
shapiro.test(butyric$value)
#W = 0.92274, p-value = 0.00099
shapiro.test(valeric$value)
#W = 0.62978, p-value = 4.883e-11
shapiro.test(caproic$value)
#W = 0.44446, p-value = 1.026e-13
ggqqplot(valeric$value) #normal tests showed none of the data are normal.

tx.total <- with(total, interaction(Sample, variable))
kruskal.test(value ~ tx.total, data=total)
# Kruskal-Wallis rank sum test
# 
# data:  value by tx.total
# Kruskal-Wallis chi-squared = 41.634, df = 9, p-value = 3.834e-06


kruskal_total<-kruskal(total$value, tx.total, group=TRUE, p.adj="BH")

tx.acetic <- with(acetic, interaction(Sample, variable))
kruskal_acetic<-kruskal(acetic$value, tx.acetic, group=TRUE, p.adj="BH")

tx.propanoic <- with(propanoic, interaction(Sample, variable))
kruskal_propanoic <- kruskal(propanoic$value, tx.propanoic, group=TRUE, p.adj="BH")

tx.butyric <- with(butyric, interaction(Sample, variable))
kruskal_butyric<-kruskal(butyric$value, tx.butyric, group=TRUE, p.adj="BH")

tx.valeric <- with(valeric, interaction(Sample, variable))
kruskal_valeric<-kruskal(valeric$value, tx.valeric, group=TRUE, p.adj="BH")

tx.caproic <- with(caproic, interaction(Sample, variable))
kruskal_caproic<-kruskal(caproic$value, tx.caproic, group=TRUE, p.adj="BH")

vfa_letters <- as.data.frame(c(kruskal_total$groups,
                               kruskal_acetic$groups, 
                               kruskal_propanoic$groups, 
                               kruskal_butyric$groups, 
                               kruskal_valeric$groups, 
                               kruskal_caproic$groups))
row.names(vfa_letters)<- row.names(kruskal_total$groups)
