setwd("path/to/working/directory")
nh3tab<-read.table("ammonia.txt", sep="\t", header=TRUE, check.names= FALSE)

library(ggplot2)
library(dplyr)
library(agricolae)
library(ggpubr)

# > sessionInfo()
# R version 3.6.0 (2019-04-26)
# [1] ggpubr_0.2.5    magrittr_1.5    agricolae_1.3-1 dplyr_0.8.1     ggplot2_3.1.1  



sample_table <- table(nh3tab$Sample)
nh3tab$Sample2 <- factor(nh3tab$Sample, levels = c("E. pyruvativorans isol6",
                                             "A. sticklandii 12662",
                                             "C. aminophilum 49906",
                                             "R. flavefaciens 007",
                                             "R. albus SY3",
                                             "F. succinogenes S85",
                                             "B. fibrisolvens 3071",
                                             "P. ruminicola 19189",
                                             "M. elsedenii T81",
                                             "M2 medium")
)

blank<- nh3tab[(grep("M2 medium", nh3tab$Sample)), ]

shapiro.test(nh3tab$Amount) #normal test showed not normal.
# Shapiro-Wilk normality test
# 
# data:  nh3tab$Amount
# W = 0.77191, p-value = 2.872e-08

ggqqplot(nh3tab$Amount) #normal test showed not normal.



tx.nh3 <- with(nh3tab, interaction(Sample))
kruskal.test(Amount ~ tx.nh3, data=nh3tab)
# Kruskal-Wallis rank sum test
# data:  Amount by tx.nh3
# Kruskal-Wallis chi-squared = 53.251, df = 9, p-value = 2.619e-08

kruskal_nh3<-kruskal(nh3tab$Amount, tx.nh3, group=TRUE, p.adj="bonferroni") #in Agricolae
letters <- kruskal_nh3$groups
bacteria_names <- unique(nh3tab$Sample)
letters <- letters[ order(match(row.names(letters), bacteria_names)), ]


nh3.summarized = nh3tab %>% group_by(Sample2) %>% summarize(Max.nh3=max(Amount))

ggplot(nh3tab, aes(x=Sample2, y=Amount, fill=nh3tab$Phenotype, label=NA)) + 
  geom_boxplot() +
  labs(x= "Organism") +
  labs(y =  expression(paste("Amount of Ammonia (  ", mu, "g/ml)"))) + #put the mu symbol in
 # ggtitle("Ammonia Production") +
  theme(text = element_text(size=16), 
        axis.text.x = element_text(angle=60, hjust=1, face="italic"),
        axis.title.y = element_text(vjust=3)
        ) +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF"), name="Phenotype") + 
  geom_hline(yintercept= median(blank$Amount)) + 
  geom_hline(yintercept= max(blank$Amount), linetype = "dashed") + 
  geom_hline(yintercept= min(blank$Amount), linetype = "dashed") +
  geom_text(inherit.aes = FALSE,
   data=nh3.summarized, 
   aes(x=nh3.summarized$Sample2, y=15+nh3.summarized$Max.nh3, label=letters$groups))
  # vjust=0.1
  #+
 # scale_y_continuous(expand = c(0, 0), limits = c(0, 450))

