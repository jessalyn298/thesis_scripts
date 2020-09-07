library(ggplot2)
library(reshape2)
library(ggforce)


#sessionInfo()
# other attached packages:
#   [1] ggforce_0.3.1  reshape2_1.4.3 ggplot2_3.1.1


#Set working directory
setwd("path/to/working/directory")

boxplot_data<-read.csv("boxplot_data.csv", sep=",", 
                       check.name=FALSE, header=TRUE)

melteddata<- melt(boxplot_data, id=c("Species", "Phenotype"))
melteddata$Species<- factor(melteddata$Species, levels = unique(melteddata$Species))

plot1<-
ggplot(melteddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=1.75, position=position_dodge2(width=0.5), colour="black") +
  geom_point(size=1.25, position=position_dodge2(width=0.5), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page = 1, nrow = 9, ncol = 5 ) +
  theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 6))



n_pages(plot1) # 6

plot2<- 
  ggplot(melteddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=1.75, position=position_dodge2(width=0.5), colour="black") +
  geom_point(size=1.25, position=position_dodge2(width=0.5), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page = 2, nrow = 9, ncol = 5 ) +
  theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 6))

plot3<- 
  ggplot(melteddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=1.75, position=position_dodge2(width=0.5), colour="black") +
  geom_point(size=1.25, position=position_dodge2(width=0.5), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page = 3, nrow = 9, ncol = 5 ) +
  theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 6))

plot4<- 
  ggplot(melteddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=1.75, position=position_dodge2(width=0.5), colour="black") +
  geom_point(size=1.25, position=position_dodge2(width=0.5), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page = 4, nrow = 9, ncol = 5 ) +
  theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 6))

plot5<- 
  ggplot(melteddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=1.75, position=position_dodge2(width=0.5), colour="black") +
  geom_point(size=1.25, position=position_dodge2(width=0.5), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page = 5, nrow = 9, ncol = 5 ) +
  theme(axis.text.x = element_blank(), strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 6))

plot6<- 
  ggplot(melteddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=1.75, position=position_dodge2(width=0.5), colour="black") +
  geom_point(size=1.25, position=position_dodge2(width=0.5), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page = 6, nrow = 9, ncol = 5 ) +
  theme(axis.text.x = element_blank())


ggsave("boxplot_1.png", device="png", plot=plot1, dpi=300, width = 16.5, height = 26, units="cm")
ggsave("boxplot_2.png", device="png", plot=plot2, dpi=300, width = 16.5, height = 26, units="cm")
ggsave("boxplot_3.png", device="png", plot=plot3, dpi=300, width = 16.5, height = 26, units="cm")
ggsave("boxplot_4.png", device="png", plot=plot4, dpi=300, width = 16.5, height = 26,units="cm")
ggsave("boxplot_5.png", device="png", plot=plot5, dpi=300, width = 16.5, height = 26, units="cm")
ggsave("boxplot_6.png", device="png", plot=plot6, dpi=300, width = 16.5, height = 26, units="cm")






####boxplots for figures###
boxplot_data_reduced<-read.csv("interesting_letters_boxplots_part2.csv", sep=",", 
                       check.name=FALSE, header=TRUE)

melted_reduceddata<- melt(boxplot_data_reduced, id=c("Species", "Phenotype"))
melted_reduceddata$Species<- factor(melted_reduceddata$Species, levels = unique(melted_reduceddata$Species))

#reduced_bp<-
plotplot1<-ggplot(melted_reduceddata, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("eggNOG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=1, nrow = 6, ncol = 4 ) +
  theme(axis.text.x = element_blank())
  
plotplot2<-ggplot(melted_reduceddata, aes(x=variable, y=value, fill=Phenotype)) +
    geom_boxplot() +
    scale_fill_manual(values=c("#F8766D", "#00BA38")) +
    geom_point(size=2, position=position_dodge2(width=1), colour="black") +
    geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
    scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
    xlab("eggNOG Functional Orthologous Group") + 
    ylab("Normalised Count Data") +
    facet_wrap_paginate(~variable, scales='free', page=2, nrow = 6, ncol = 4) +
    theme(axis.text.x = element_blank())
    #theme(strip.background = element_blank(), strip.text = element_blank())

n_pages(reduced_bp) #2

ggsave("figures_boxplot_1.png", device="png", plot=plotplot1, dpi=300, width = 18, height = 25, units="cm")
ggsave("figures_boxplot_2.png", device="png", plot=plotplot2, dpi=300, width = 18, height = 25, units="cm")