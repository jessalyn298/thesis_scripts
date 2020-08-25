library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggforce)

#Set working directory
setwd("path/to/working/directory")

#get data
boxplotdat<-read.table("interesting_letters_for_boxplot.csv", sep=",", header=TRUE, check.names=FALSE)

#make data look like this:
#cluster species phenotype count
boxplotdat_melted<-melt(boxplotdat, id=c("Species", "Phenotype"))
#force ordering so x-axis labels don't get moved around
boxplotdat_melted$Species<- factor(boxplotdat_melted$Species, levels = unique(boxplotdat_melted$Species))


###make lots of boxplots to see which ones are interesting###
# bp1<-
#   ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
#   geom_boxplot() +
#   scale_fill_manual(values=c("#F8766D", "#00BA38")) +
#   geom_point(size=2, position=position_dodge2(width=1), colour="black") +
#   geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
#   scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
#   xlab("KEGG Functional Orthologous Group") + 
#   ylab("Normalised Count Data") +
#   facet_wrap_paginate(~variable, scales='free', page=1, nrow = 5, ncol = 6 ) +
#   theme(axis.text.x = element_blank())
#n_pages(bp1) # 9

plot1<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=1, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot2<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=2, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot3<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=3, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot4<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=4, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot5<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=5, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot6<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=6, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot7<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=7, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot8<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=8, nrow = 5, ncol = 6 ) +
  theme(axis.text.x = element_blank())

plot9<- 
  ggplot(boxplotdat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=9, nrow = 3, ncol = 6 ) +
  theme(axis.text.x = element_blank())

# plot1
# plot2
# plot3
# plot4
# plot5
# plot6
# plot7
# plot8
# plot9

#save them as png images:

ggsave("boxplot_1.png", device="png", plot=plot1, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_2.png", device="png", plot=plot2, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_3.png", device="png", plot=plot3, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_4.png", device="png", plot=plot4, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_5.png", device="png", plot=plot5, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_6.png", device="png", plot=plot6, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_7.png", device="png", plot=plot7, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_8.png", device="png", plot=plot8, dpi=300, width = 25, height = 18, units="cm")
ggsave("boxplot_9.png", device="png", plot=plot9, dpi=300, width = 25, height = 18, units="cm")




###Make boxplots for those interesting ones:

interesting_dat<-read.table("interesting_letters_for_boxplot_reduced.csv", sep=",", header=TRUE, check.names=FALSE)
interesting_dat_melted<-melt(interesting_dat, id=c("Phenotype", "Species"))
interesting_dat_melted$Species<- factor(interesting_dat_melted$Species, levels = unique(interesting_dat_melted$Species))

#reduced_bp<-
plotplot1<-ggplot(interesting_dat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=1, nrow = 6, ncol = 5) +
  theme(axis.text.x = element_blank())

plotplot2<-ggplot(interesting_dat_melted, aes(x=variable, y=value, fill=Phenotype)) +
  geom_boxplot() +
  scale_fill_manual(values=c("#F8766D", "#00BA38")) +
  geom_point(size=2, position=position_dodge2(width=1), colour="black") +
  geom_point(size=1.3, position=position_dodge2(width=1), aes(colour=Species)) +
  scale_color_manual(values=c("#FF8000", "#FF0000", "#FF80C0", "#80FF00", "#008000", "#3A9494")) +
  xlab("KEGG Functional Orthologous Group") + 
  ylab("Normalised Count Data") +
  facet_wrap_paginate(~variable, scales='free', page=2, nrow = 6, ncol = 5) +
  theme(axis.text.x = element_blank())
#theme(strip.background = element_blank(), strip.text = element_blank())


ggsave("K0s_figures_boxplot_1.png", device="png", plot=plotplot1, dpi=300, width = 18, height = 25, units="cm")
ggsave("K0s_figures_boxplot_2.png", device="png", plot=plotplot2, dpi=300, width = 18, height = 25, units="cm")