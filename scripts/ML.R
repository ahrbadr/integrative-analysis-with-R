
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# # Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

library("ggplot2")
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(dplyr)
library(readr)






##### load the mRNA-Seq data #####

setwd("E:/Integrative/smoc2")


smoc2_rawcounts <- read.csv("fibrosis_smoc2_rawcounts_unordered.csv",row.names = 1)

# Explore the first five observations of smoc2_rawcounts
head(smoc2_rawcounts)

# Explore the structure of smoc2_rawcounts
str(smoc2_rawcounts)

#Use the information below to create a metadata data frame for the fibrosis count data called metadata with columns genotype and condition. The sample names (e.g. smoc2_fibrosis1, smoc2_fibrosis2, etc.) should be the row names of the data frame:

###############
# The MetData #
###############
# Create genotype vector
genotype <- c("smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe", "smoc2_oe")

# Create condition vector
condition <- c("fibrosis", "fibrosis", "fibrosis", "fibrosis", "normal", "normal", "normal")

# Create data frame
smoc2_metadata <- data.frame(genotype, condition)

# Assign the row names of the data frame
rownames(smoc2_metadata) <- c("smoc2_fibrosis1", "smoc2_fibrosis2", "smoc2_fibrosis3", "smoc2_fibrosis4", "smoc2_normal1", "smoc2_normal3", "smoc2_normal4")

# Matching metadata and counts data


# Use the match() function to reorder the columns of the raw counts
reorder_idx <- match(rownames(smoc2_metadata), colnames(smoc2_rawcounts))

# Reorder the columns of the count data
reordered_smoc2_rawcounts <- smoc2_rawcounts[ , reorder_idx]

# Create a DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                                    colData = smoc2_metadata,
                                    design = ~ condition)

# 
# Normalizing counts with DESeq2

# Determine the size factors to use for normalization
dds_smoc2 <- estimateSizeFactors(dds_smoc2)
sizeFactors(dds_smoc2)

# Extract the normalized counts
smoc2_normalized_counts <- counts(dds_smoc2, normalized=TRUE)

#####################################
#Hierarchical heatmap by condition
#######################################


#When performing quality assessment for count data,
#we need to transform the normalized counts for better visualization of the variance for unsupervised clustering analyses. 
#To assess the similarity of the smoc2 samples using hierarchical heatmaps, transform the normalized counts and perform hierarchical clustering analysis. 

# Transform the normalized counts 
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)

# Compute the correlation values between samples
vsd_cor_smoc2 <- cor(vsd_mat_smoc2) 
#View correlation values
View(vsd_cor_smoc2) 

# The annotation argument selects which factors in the
# metadata to include as annotation bars. We use the select()
# function from the dplyr package to select the condition 
# column in the smoctype metadata

# Plot the heatmap
pheatmap(vsd_cor_smoc2, annotation = select(smoc2_metadata, condition))

#The biological replicates cluster together 
# and the samples in different conditions cluster separately.
# There are no outliers or samples with low correlation values relative to all other samples.






#########################################

#########################################
# PCA analysis
# To continue with the quality assessment of our samples, in the first part of this exercise, we will perform PCA to look how our samples cluster and whether our condition of interest corresponds with the principal components explaining the most variation in the data. In the second part, we will answer questions about the PCA plot.
# 
# To assess the similarity of the smoc2 samples using PCA, we need to transform the normalized counts then perform the PCA analysis. Assume all libraries have been loaded, the DESeq2 object created, and the size factors have been stored in the DESeq2 object, dds_smoc2.
# 
# Run the code to transform the normalized counts.
# Perform PCA by plotting PC1 vs PC2 using the DESeq2 plotPCA() function on the DESeq2 transformed counts object, vsd_smoc2 and specify the intgroup argument as the factor to color the plot.
# 

# Transform the normalized counts 
vsd_smoc2 <- vst(dds_smoc2, blind = TRUE)

# Plot the PCA of PC1 and PC2
plotPCA(vsd_smoc2, intgroup="condition")




#Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components.
#This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.
pcaData <- plotPCA(vsd_smoc2, intgroup=c("condition", "genotype"), returnData=TRUE)
percentVar <- round( 100 *attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()



# Extract the matrix of transformed counts
vsd_mat_smoc2 <- assay(vsd_smoc2)
pc.cr <- princomp(vsd_mat_smoc2, cor = TRUE)
screeplot(pc.cr)


#Related to the distance matrix is the PCA plot, which shows the samples in the 2D plane spanned by their first two principal components. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.
pcaData <- plotPCA(vsd_smoc2, intgroup=c("condition", "genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


# Comment on PCA: 
#   Real datasets come with missing values and unfortunately PCA cannot
# deal with missing values and observations containing NA will be dropped automatically.Additionally PCA
# can also not handle categorical data, 
#therefore these can be encoded as binary dummy variables. such that
# there is a distance between the categories.
# Create DESeq2 object
dds_smoc2 <- DESeqDataSetFromMatrix(countData = reordered_smoc2_rawcounts,
                                    colData = smoc2_metadata,
                                    design = ~ condition)

# Run the DESeq2 analysis
dds_smoc2 <- DESeq(dds_smoc2)


# DESeq2 model - exploring dispersions

# Plot dispersions
plotDispEsts(dds_smoc2)

#   Good fit, since the dispersions decrease with increasing mean and cluster around the maximum likelihood (ML) line.
# The assumptions of DESeq2 are met since the dispersions decrease with increasing mean and the raw dispersions seem to cluster around the maximum likelihood line.


# DESeq2 model - extracting results

# Extract the results of the differential expression analysis
smoc2_res <- results(dds_smoc2, 
                     contrast = c("condition", "fibrosis", "normal"), 
                     alpha = 0.05)

# Extract results
smoc2_res <- results(dds_smoc2, 
                     contrast = c("condition", "fibrosis", "normal"), 
                     alpha = 0.05, 
                     lfcThreshold = 0.32)

# Shrink the log2 fold changes
smoc2_res <- lfcShrink(dds_smoc2, 
                       contrast = c("condition", "fibrosis", "normal"), 
                       res = smoc2_res)

##########################################



# remove nulls
smoc2_res=smoc2_res[complete.cases(smoc2_res), ]
summary(smoc2_res)


res.df=as.data.frame(smoc2_res)
res.degs=res.df[res.df$padj< 0.01 & abs(res.df$log2FoldChange)>log2(4),]

#expression of these degs for plots

exp.degs= smoc2_rawcounts[rownames(smoc2_rawcounts) %in% rownames(res.degs), ]


# Subset normalized counts to significant genes
sig_norm_counts_smoc2 <- smoc2_normalized_counts[rownames(smoc2_res_sig), ]

# Choose heatmap color palette
heat_colors <- brewer.pal(n = 6, name = "YlOrRd")

# The annotation argument selects which factors in the
# metadata to include as annotation bars. We use the select()
# function from the dplyr package to select the condition 
# column in the smoctype metadata
# Plot heatmap
pheatmap(sig_norm_counts_smoc2, 
         color = heat_colors, 
         cluster_rows = TRUE, 
         show_rownames = FALSE,
         annotation = select(smoc2_metadata, condition), 
         scale = "row")
#########################
# Unsupervised Learning #
#########################


#Hierachical clustering


#Clustering for Samples
d <- dist(t(exp.degs))
hcl <- hclust(d)
hcl
plot(hcl)
abline(h = 5e+05, col = "red")

# 
# When considering the dendrogram, the number of clusters is defined by the distance. 
#At a specific decided
# distance the dendrogram is cut and the number of cluster corresponds to the number of remaining leafes in
# the dendrogram (tree).



sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)





# Plot a line to show the cut
abline(h = 15, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]



###################################################################

# k-means clusteirng
# The k- means clusterin algorithm aims at partitioning n observations into a fixed number of k clusters. In R
# the k- means algorithm is already provided in the stats package.
# ?kmeans
# The parameters of k-means are :
#   . X numeric data matrix
# . centers are the pre-defined number of clusters
# . nstart gives the number of times, 
#the k-means algorithm is repeated, because it has a random component.
# By repetition the model improves.





?kmeans
par(mfrow = c(1, 1))

# subset normal1 vs fibrosis1
i <- grep("1", colnames(vsd_mat_smoc2))
x <- vsd_mat_smoc2[, i]
cl <- kmeans(x, 2, nstart = 10)
plot(x, col = cl$cluster)

## With all data


x<-vsd_mat_smoc2
cl <- kmeans(x, 3, nstart = 10)
plot(x, col = cl$cluster)


# subset fibrosis

e <- grep("fibrosis", names(exp.degs))
y <- exp.degs[, e]
cl <- kmeans(y, 3, nstart = 10)
plot(y, col = cl$cluster)

# subset normal

e <- grep("normal", names(exp.degs))
y <- exp.degs[, e]
cl <- kmeans(y, 3, nstart = 10)
plot(y, col = cl$cluster)




