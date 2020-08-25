library(ggplot2)
library(reshape2)
library(ggsignif)


setwd("path/to/working/directory")
dat <- read.delim("aa_transporter_counts.tsv", check.names=FALSE)


dat1<- melt(dat, id=c("species", "phenotype"))
#This 'melts' the data, turning it into a table like this:
#bacteria 1 group 1 gene 1 3
#bacteria 1 group 1 gene 2 1
#bacteria 2 group 1 gene 1 1
#bacteria 2 group 1 gene 2 5
#etc

#to assign colours (I need 21 colours as I have 21 transporters)
mycols <-c("navy", "aquamarine", "cadetblue3", "deepskyblue4", "darkorchid4",
           "coral", "burlywood", "chocolate",  
           "firebrick1", "lightpink", "maroon4", "deeppink2", 
           "beige",
           "darkolivegreen3", "limegreen", "seagreen",
           "black") 

ggplot(dat1, aes(x = species, y = value, fill = variable)) +
  geom_bar(stat = "identity")+
  facet_wrap(~phenotype, scales = "free_x") +
  scale_fill_manual(values = mycols)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, face="italic"),
        plot.margin = margin(10, 10, 10, 55),
        axis.title.y=element_text(margin=margin(t = 0, r = 10, b = 0, l = 0)),
        text = element_text(size=12)) +
  labs(x = "Bacteria") +
  labs(y = "Number of Transporter Family Orthologs") +
  #ggtitle("Number of Amino Acid Transporter Gene Orthologs for Three Different Groups")+
  guides(fill=guide_legend(title="Transporter Families")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50))

ggsave("aatransporters_stackedgraph.png", device="png", dpi = 300)

# x = as.factor(total$Group)
# 
# ggplot(total, aes(x=Group, y=Total, fill=Group))+
#   geom_boxplot() +
#   geom_dotplot(aes(fill=Group),binaxis="y",stackdir="center",dotsize=0.4) +
#   geom_signif(comparisons = list(c("HAP", "NAP")), y_position=c(29),
#               map_signif_level=TRUE) +
#   geom_signif(comparisons = list(c("SAP", "NAP")), y_position=c(31),
#               map_signif_level=TRUE) +
#   geom_signif(comparisons = list(c("SAP", "HAP")), y_position=c(33),
#               map_signif_level=TRUE) +
#   scale_color_manual(values=c("red", "green", "blue")) +
#   labs(x = "Phenotypic Group") +
#   labs(y = "Total Number of Amino Acid and Peptide Transporter Families") +
#   ggtitle("Amino Acids and Peptide Transporter Families")+
#   guides(fill=FALSE)
# 
# 
# HAP <- total[grep ("HAP", total$Group),]
# SAP <- total[grep ("SAP", total$Group),]
# NAP <- total[grep ("NAP", total$Group),]
# 
# wilcox.test(HAP$Total, SAP$Total)
# wilcox.test(SAP$Total, NAP$Total)
# wilcox.test(HAP$Total, NAP$Total)