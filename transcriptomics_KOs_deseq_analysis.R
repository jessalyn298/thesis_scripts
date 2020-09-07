library(DESeq2)
library(pheatmap)
library(igraph)


#Set working directory
setwd("path/to/working/directory")

#get data
count_data<-read.table("../K0_counts_RNA_data_reduced.csv", header=TRUE, sep=",", row.names=1, check.names=FALSE)

#set the species sampled
sampledSpecies<-c( "12662", "49906", "Isol6", "S85", "007c", "Sy3")

#apply this to each of the columns in the data - it must match
AllSpecies <- c( "12662", "12662", "12662", "12662", "49906", "49906", "49906", "49906", "Isol6", "Isol6", "Isol6", "Isol6", "S85", "S85", "S85", "S85", "007c", "007c", "007c", "007c", "Sy3", "Sy3", "Sy3", "Sy3")

#set the groups that the species belong to
sampledGroups<-c("HAP", "NAP")

#apply this to each of the columns in the data- again it must match
AllGroup <- c("HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "HAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP", "NAP")

#set lower limit for the total sum of the rows
#gt1<- function(x){
#  sum(x) > 0;
#}

#Set definitions
#MIN_SAMPLE_READS = 1000000
#MIN_GENE_READS = 100

#only include samples that have more than the minimum number of reads specified above
#aboveminsize<-count_data[,apply(count_data, 2, sum) > MIN_SAMPLE_READS]

#Only include genes that have more than the minimum number of reads defined above.#
#aboveminreads = aboveminsize[apply((aboveminsize), 1, sum) > MIN_GENE_READS,]

# get rid of all genes that present a zero in any sample
#reduced_count_data=aboveminreads[apply(aboveminreads, 1, min) > 0,]

#data_for_analysis = aboveminreads # set this to which of the datasets is to be analysed.

# Do DESEQ2 Bit
coldata = data.frame(
  Species= AllSpecies,
  row.names=colnames(count_data)
)

dds <- DESeqDataSetFromMatrix(
  countData = count_data,
  colData = coldata,
  design = ~ Species
)

ddsMF <- DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# -- note: fitType='parametric', but the dispersion trend was not well captured by the
# function: y = a/x + b, and a local regression fit was automatically substituted.
# specify fitType='local' or 'mean' to avoid this message next time.
# final dispersion estimates
# fitting model and testing

resultsNames(ddsMF)
# [1] "Intercept"             "Species_12662_vs_007c" "Species_49906_vs_007c" "Species_Isol6_vs_007c"
# [5] "Species_S85_vs_007c"   "Species_Sy3_vs_007c"  


# Generate normalsied values	
dds <- estimateSizeFactors(dds)
normalised_counts<-counts(dds, normalized=TRUE)

write.table(normalised_counts, "normalised_counts_rna_data_K0s.table", sep="\t")

# Generate statistical results for all-against-all comparisons of species.
# Species to compare: "12662", "49906", "Isol6", "S85", "007c", "Sy3"

Res12662_49906 <- results(ddsMF, contrast=c("Species", "12662", "49906" ))
sum(Res12662_49906$padj < 0.1, na.rm=TRUE) #507

Res12662_isol6 <- results(ddsMF, contrast=c("Species", "12662", "Isol6" ))
sum(Res12662_isol6$padj < 0.1, na.rm=TRUE) #511

Res12662_S85 <- results(ddsMF, contrast=c("Species", "12662", "S85" ))
sum(Res12662_S85$padj < 0.1, na.rm=TRUE) #515

Res12662_007c <- results(ddsMF, contrast=c("Species", "12662", "007c" ))
sum(Res12662_007c$padj < 0.1, na.rm=TRUE) #504

Res12662_Sy3 <- results(ddsMF, contrast=c("Species", "12662", "Sy3" ))
sum(Res12662_Sy3$padj < 0.1, na.rm=TRUE) #534

Res49906_Isol6 <- results(ddsMF, contrast=c("Species", "49906", "Isol6" ))
sum(Res49906_Isol6$padj < 0.1, na.rm=TRUE) #500

Res49906_S85 <- results(ddsMF, contrast=c("Species", "49906", "S85" ))
sum(Res49906_S85$padj < 0.1, na.rm=TRUE) #517

Res49906_007c <- results(ddsMF, contrast=c("Species", "49906", "007c" ))
sum(Res49906_007c$padj < 0.1, na.rm=TRUE) #504

Res49906_Sy3 <- results(ddsMF, contrast=c("Species", "49906", "Sy3" ))
sum(Res49906_Sy3$padj < 0.1, na.rm=TRUE) #535

ResIsol6_S85 <- results(ddsMF, contrast=c("Species", "Isol6", "S85" ))
sum(ResIsol6_S85$padj < 0.1, na.rm=TRUE) #503

ResIsol6_007c <- results(ddsMF, contrast=c("Species", "Isol6", "007c" ))
sum(ResIsol6_007c$padj < 0.1, na.rm=TRUE) #511

ResIsol6_Sy3 <- results(ddsMF, contrast=c("Species", "Isol6", "Sy3" ))
sum(ResIsol6_Sy3$padj < 0.1, na.rm=TRUE) #523

ResS85_007c <- results(ddsMF, contrast=c("Species", "S85", "007c" ))
sum(ResS85_007c$padj < 0.1, na.rm=TRUE) #508

ResS85_Sy3 <- results(ddsMF, contrast=c("Species", "S85", "Sy3" ))
sum(ResS85_Sy3$padj < 0.1, na.rm=TRUE) #520

Res007c_Sy3 <- results(ddsMF, contrast=c("Species", "007c", "Sy3" ))
sum(Res007c_Sy3$padj < 0.1, na.rm=TRUE) #523


#compile the padjusted significant results into one table
overall_results<-cbind(Res12662_49906[,6], Res12662_isol6[,6], Res12662_S85[,6], Res12662_Sy3[,6], Res12662_007c[,6], Res49906_Isol6[,6], Res49906_S85[,6], Res49906_Sy3[,6], Res49906_007c[,6], ResIsol6_S85[,6], ResIsol6_Sy3[,6], ResIsol6_007c[,6], ResS85_Sy3[,6], ResS85_007c[,6], Res007c_Sy3[,6])
colnames(overall_results)<-c("12662_49906", "12662_Isol6", "12662_S85", "12662_Sy3", "12662_007c", "49906_Isol6", "49906_S85", "49906_Sy3", "49906_007c", "Isol6_S85", "Isol6_Sy3", "ISol6_007c", "S85_Sy3", "S85_007c", "007c_Sy3")
rownames(overall_results)<-rownames(count_data)

#Compile the logfold results into one table
overall_logfold<-cbind(Res12662_49906[,2], Res12662_isol6[,2], Res12662_S85[,2], Res12662_Sy3[,2], Res12662_007c[,2], Res49906_Isol6[,2], Res49906_S85[,2], Res49906_Sy3[,2], Res49906_007c[,2], ResIsol6_S85[,2], ResIsol6_Sy3[,2], ResIsol6_007c[,2], ResS85_Sy3[,2], ResS85_007c[,2], Res007c_Sy3[,2])
colnames(overall_logfold)<-c("12662_49906", "12662_Isol6", "12662_S85", "12662_Sy3", "12662_007c", "49906_Isol6", "49906_S85", "49906_Sy3", "49906_007c", "Isol6_S85", "Isol6_Sy3", "ISol6_007c", "S85_Sy3", "S85_007c", "007c_Sy3")
rownames(overall_logfold)<-rownames(count_data)

#set all NAs to 1
overall_results[is.na(overall_results)] <- 1

write.table(overall_results, "overall_results_HapNapRNA_allcomparisons_K0s.txt", sep='\t')

# This function takes a single line of the overall_results talbe at a time and calculates the signifincance groups and outputs those in letter format
generate_significance_groups<-function(pvals)
{  
  #Function to calculate significance groups, pass this a line from the overall_results table
  graph<-NULL #define an empty graph

  # for each of the consolidated significance results generate an edge in the graph between Time points which are NOT significantly different from each other
  if(pvals[1] >= 0.1) graph<-c(graph, "12662", "49906")
  if(pvals[2] >= 0.1) graph<-c(graph, "12662", "Isol6")
  if(pvals[3] >= 0.1) graph<-c(graph, "12662", "S85")
  if(pvals[4] >= 0.1) graph<-c(graph, "12662", "Sy3")
  if(pvals[5] >= 0.1) graph<-c(graph, "12662", "007c")
  if(pvals[6] >= 0.1) graph<-c(graph, "49906", "Isol6")
  if(pvals[7] >= 0.1) graph<-c(graph, "49906", "S85")
  if(pvals[8] >= 0.1) graph<-c(graph, "49906", "Sy3")
  if(pvals[9] >= 0.1) graph<-c(graph, "49906", "007c")
  if(pvals[10] >= 0.1) graph<-c(graph, "Isol6", "S85")	
  if(pvals[11] >= 0.1) graph<-c(graph, "Isol6", "Sy3")
  if(pvals[12] >= 0.1) graph<-c(graph, "Isol6", "007c")
  if(pvals[13] >= 0.1) graph<-c(graph, "S85", "Sy3")
  if(pvals[14] >= 0.1) graph<-c(graph, "S85", "007c")
  if(pvals[15] >= 0.1) graph<-c(graph, "007c", "Sy3")

  # Generate self-loops in the graph to ensure that all timepoints are represented in the graph
  graph<-c(graph, "12662", "12662", "49906", "49906", "Isol6", "Isol6", "S85", "S85","Sy3", "Sy3", "007c", "007c" )
  
  # Generate the graph from the edges defined in g
  g<-graph(graph)
  #plot(g)

  #Find the cliques in the graph to identify the significance groups
  r<-suppressWarnings(max_cliques(g))
  
  #Turn signifncance groups into a matrix representation for allow printing of signifncance groups with letters.
  #This will contain rowsa for each clique (signifiance group), and colums for each sample, with '1' where a sample is a member of a clique and '0' where it is not.
  a<-matrix(0,nrow=length(r),ncol=length(sampledSpecies));
  colnames(a)<-sampledSpecies
  rownames(a)<-letters[1:length(r)]
  a[,rownames(as.matrix(unlist(r[1])))][1]<-1
  for(j in 1:length(r)){
    for(i in rownames(as.matrix(unlist(r[j]))) ) {
      a[,i][j]<-1
    }
  }
  
  #This generates a vector which contains letters for each sample, detailing their membership in signifcance groups.
  consolidate<-matrix(NA, nrow=1, ncol=length(sampledSpecies))
  colnames(consolidate)<-sampledSpecies
  for (i in 1:length(sampledSpecies)){
    consolidate[1, sampledSpecies[i]]=paste(rownames(a)[grep("1", a[, sampledSpecies[i]])], collapse=" ")
  }
  
  # The for loop above automates the following lines
  #consolidate[1,'T1']=paste(rownames(a)[grep("1", a[,'T1'])], collapse=" ")
  #consolidate[1,'T2']=paste(rownames(a)[grep("1", a[,'T2'])], collapse=" ")
  #consolidate[1,'T4']=paste(rownames(a)[grep("1", a[,'T4'])], collapse=" ")
  #consolidate[1,'T6']=paste(rownames(a)[grep("1", a[,'T6'])], collapse=" ")
  #consolidate[1,'T8']=paste(rownames(a)[grep("1", a[,'T8'])], collapse=" ")
  return(consolidate)
}

# Only include genes that were significant for the "post hoc" tests
#Send all results to file (only significant)
overall_results_sig_only<- overall_results [apply(overall_results, 1, min) < 0.1,]
write.table(overall_results_sig_only, "Overall_results_sig_only_HapNapRNA_k0.txt", sep='\t')

#Send logfold results to file
overall_logfold_sig_only<-overall_logfold[apply(overall_results, 1, min) < 0.1,]
write.table(overall_logfold_sig_only, "Overall_logFold_sig_only_HAPNAPRNA_K0_txt", sep='\t')

# Generate the matrix of normalised counts for those genes that are retained for significance
HapNapRNA_normalised_counts_Sigonly<-normalised_counts[rownames(overall_results_sig_only),]
  
# Calculate mean values for each group or treatment
mean_norm_counts<-matrix(NA, ncol=length(sampledSpecies), nrow=dim(HapNapRNA_normalised_counts_Sigonly)[1])
colnames(mean_norm_counts)<-sampledSpecies
rownames(mean_norm_counts)<-rownames(HapNapRNA_normalised_counts_Sigonly)
for(i in 1:length(sampledSpecies))
{
  mean_norm_counts[,sampledSpecies[i]]<-apply(HapNapRNA_normalised_counts_Sigonly[,grep(sampledSpecies[i], AllSpecies)], 1, mean)
}

write.table(mean_norm_counts, "mean_norm_counts_HapNapRNA_sig_genes_K0.txt", sep='\t')

# Create the matrix to hold the results of the significance groups analysis
sig_groups<-matrix(NA, nrow=dim(HapNapRNA_normalised_counts_Sigonly)[1], ncol=length(sampledSpecies))
colnames(sig_groups)<-sampledSpecies
rownames(sig_groups)<-rownames(HapNapRNA_normalised_counts_Sigonly)

# Calculate the signifncance groups for every sig gene
for(i in 1:dim(overall_results_sig_only)[1])
{
  sig_groups[i,]<-generate_significance_groups(overall_results_sig_only[i,])
}

write.table(sig_groups, "K0_significance_letters.txt", sep='\t')

  
  
