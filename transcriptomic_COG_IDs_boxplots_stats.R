library(ggplot2)
library(reshape2)
library(ggpubr)

#[1] ggpubr_0.2.5   reshape2_1.4.3 ggplot2_3.1.1 

setwd("path/to/working/directory")
#dat <- read.delim("functional_category_percentage_counts.txt", check.names=FALSE)
dat<-read.delim("functional_category_percentage_counts_bp_data.txt", check.names=FALSE)

#I don't need to do normality tests as I know that my data isn't normal - it's proportional and bounded (0-100%).
list_of_categories <- colnames(dat)[-c(1,2)]
for (i in list_of_categories) {
  print(i)
  print(shapiro.test(as.matrix(bp_data[i])))
}

dat_HAPs <- dat[grep("HAP", dat$phenotype),]
dat_NAPs <- dat[grep("NAP", dat$phenotype),]
# 
# for (i in list_of_categories) {
#   print(i)
#   print(sd(as.matrix(dat_HAPs[i])))
#   print(sd(as.matrix(dat_NAPs[i])))
# }


for (i in list_of_categories) {
  print(i)
  print(wilcox.test(as.matrix(dat_HAPs[i]),as.matrix(dat_NAPs[i]), paired = FALSE, alternative = "two.sided" ))
}

# [1] "J"
#
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 22, p-value = 0.002914
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "A"
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 48, p-value = 0.03669
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "K"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 43, p-value = 0.1005
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "L"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 55, p-value = 0.3474
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "B"
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 72, p-value = NA
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "D"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 13, p-value = 0.0002744
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "Y"
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 72, p-value = NA
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "V"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 20, p-value = 0.00183
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "T"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 54, p-value = 0.3186
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "M"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 48, p-value = 0.1782
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "N"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 62, p-value = 0.5899
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "Z"
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 72, p-value = NA
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "W"
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 72, p-value = NA
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "U"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 0, p-value = 7.396e-07
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "O"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 29, p-value = 0.01209
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "C"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 144, p-value = 7.396e-07
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "G"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 25, p-value = 0.00556
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "E"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 126, p-value = 0.001115
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "F"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 95, p-value = 0.1978
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "H"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 77, p-value = 0.7987
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "I"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 128, p-value = 0.000656
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "P"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 141, p-value = 5.177e-06
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "Q"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 74, p-value = 0.9323
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "R"
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 72, p-value = NA
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "S"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 96, p-value = 0.1782
# alternative hypothesis: true location shift is not equal to 0
# 
# [1] "no_identifier"
# 
# Wilcoxon rank sum test
# 
# data:  as.matrix(dat_HAPs[i]) and as.matrix(dat_NAPs[i])
# W = 0, p-value = 7.396e-07
# alternative hypothesis: true location shift is not equal to 0


bp_data<-read.delim("functional_category_percentage_counts_bp_data.txt", check.names=FALSE)
bp_dat_melt<-melt(bp_data, id=c("Species", "phenotype"))
x <- as.factor(dat$Species)

dummy <- read.delim("maximum_dummydata.txt", check.names=FALSE)
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
