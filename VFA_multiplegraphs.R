setwd("path/to/working/directory")
tab<-read.table("VFA_means_sd_mM.txt", sep="\t", header=TRUE, check.names=FALSE)
#input is  table containing means and standard deviation for each species and each VFA and total.

library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

# > sessionInfo()
# R version 3.6.0 (2019-04-26)
#   [1] scales_1.0.0  gridExtra_2.3 ggplot2_3.1.1


#to get the colours of the default ggplot red, green and blue, use this command: 
show_col(hue_pal()(3))

tab$sample <- factor(tab$sample, levels = tab$sample) #THIS STOPS X AXIS REORDERING

acetate <- ggplot(tab, aes(x=sample, y=acetateMEAN, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin=acetateMEAN - acetateSD, ymax=acetateMEAN + acetateSD), width=.1) +
  labs(x= NULL) +
  labs(y = NULL) +
  ggtitle("Acetate") +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF"))+
  theme(text = element_text(size=12), 
        axis.text.x = element_text(angle=60, hjust=1, face="italic"), 
        plot.margin=margin(10,10,10,50)) + 
  geom_hline(yintercept= 20.03) + 
  geom_hline(yintercept= 25.03, linetype = "dashed") + 
  geom_hline(yintercept= 15.03, linetype = "dashed") +
  theme(legend.position="none")


propionate <- ggplot(tab, aes(x=sample, y=propionateMEAN, fill=group)) + 
  geom_bar(stat="identity", colour="black") +
  geom_errorbar(aes(ymin= tab$propionateMEAN - tab$propionateSD, ymax=propionateMEAN + propionateSD), width=.1) +
  labs(x= NULL) +
  labs(y = NULL) +
  ggtitle("Proprionate") +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF")) +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=60, hjust=1, face="italic"), plot.margin=margin(10,10,10,50)) + 
  geom_hline(yintercept= 5.13) + 
  geom_hline(yintercept= 5.64, linetype = "dashed") + 
  geom_hline(yintercept= 4.62, linetype = "dashed") +
  theme(legend.position="none")

  

butyrate <- ggplot(tab, aes(x=sample, y=butyrateMEAN, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin= butyrateMEAN - butyrateSD, ymax=butyrateMEAN + butyrateSD), width=.1) +
  labs(x= NULL) +
  labs(y = NULL) +
  ggtitle("Butyrate") +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF")) +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=60, hjust=1, face="italic"), plot.margin=margin(10,10,10,50)) +
  geom_hline(yintercept= 10.47) + 
  geom_hline(yintercept= 12.54, linetype = "dashed") + 
  geom_hline(yintercept= 8.4, linetype = "dashed") +
  theme(legend.position="none")


valerate <- ggplot(tab, aes(x=sample, y=valerateMEAN, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin= valerateMEAN - valerateSD, ymax=valerateMEAN + valerateSD), width=.1) +
  labs(x= NULL) +
  labs(y = NULL) +
  ggtitle("Valerate") +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF")) +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=60, hjust=1, face="italic"), plot.margin=margin(10,10,10,50)) + 
  geom_hline(yintercept= 3.39) + 
  geom_hline(yintercept= 3.87, linetype = "dashed") + 
  geom_hline(yintercept= 2.91, linetype = "dashed") +
  theme(legend.position="none")


caproate <- ggplot(tab, aes(x=sample, y=caproateMEAN, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin= caproateMEAN - caproateSD, ymax=caproateMEAN + caproateSD), width=.1) +
  labs(x= NULL) +
  labs(y = NULL) +
  ggtitle("Caproate") +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF")) +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=60, hjust=1, face="italic"), plot.margin=margin(10,10,10,50)) +
  geom_hline(yintercept= 4.85) + 
  geom_hline(yintercept= 5.79, linetype = "dashed") + 
  geom_hline(yintercept= 3.91, linetype = "dashed") +
  theme(legend.position="none") 


total <- ggplot(tab, aes(x=sample, y=totalMEAN, fill=group)) + 
  geom_bar(stat="identity", color="black") +
  geom_errorbar(aes(ymin= totalMEAN - totalSD, ymax=totalMEAN + totalSD), width=.1) +
  labs(x= NULL) +
  labs(y = NULL) +
  ggtitle("Total VFA Production") +
  scale_fill_manual(values=c("grey", "#F8766D", "#00BA38", "#619CFF"), name="Phenotype") +
  theme(text = element_text(size=12), axis.text.x = element_text(angle=60, hjust=1, face="italic"), legend.position="bottom", plot.margin=margin(10,10,10,50)) +
  geom_hline(yintercept= 43.86) + 
  geom_hline(yintercept= 52.57, linetype = "dashed") + 
  geom_hline(yintercept= 35.15, linetype = "dashed") 
 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

my_legend <- g_legend(total)


grid.arrange(arrangeGrob(acetate + theme(legend.position="none"),
                         propionate + theme(legend.position="none"),
                         butyrate + theme(legend.position="none"),
                         valerate + theme(legend.position="none"),
                         caproate + theme(legend.position="none"),
                         total + theme(legend.position="none"),
                         ncol=2),
             my_legend, nrow=2,heights=c(10, 1),
             bottom=textGrob("Organisms", gp=gpar(fontsize=14)),
             left=textGrob("Amount (mM)", gp=gpar(fontsize=14), rot = 90)
)
