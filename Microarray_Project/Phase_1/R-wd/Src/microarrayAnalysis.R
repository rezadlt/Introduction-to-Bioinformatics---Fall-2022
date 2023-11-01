setwd("C:/Users/Narges.Narges/Desktop/Mohammad/BIO/Project/phase1/R-wd/")

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(Rtsne)


series <- "GSE48558"
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
if (length(gset) > 1) idx <- grep("GPL6244", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#### all meta data like Group,Accession,Title,Source name,Phenotype are stored in samples_meta_table.csv
meta_table <- read.csv("samples_meta_table.csv")
gr <- c(rep("aml", 40), "healthy", rep("aml", 3), "healthy", rep("aml", 23), "healthy", "aml", "healthy", rep("aml", 3), "healthy", "aml", rep("healthy", 4), "aml", "healthy", rep("aml", 2), rep("healthy", 2), rep("aml", 2), rep("healthy", 2), "aml", "healthy", "aml", "healthy", "aml", "healthy", "aml", "healthy", "aml", "healthy", rep("aml", 3), "healthy", rep("aml", 3), "healthy", rep("aml", 29), rep("healthy", 7), rep("aml", 2), "healthy", rep("aml", 3), rep("healthy", 20))
ex <- exprs(gset)


#### drawing box plot for quality control
pdf("Results/boxplot.pdf",width = 64)
boxplot(ex)
dev.off()


#### drawing correlation heat map for all columns of data
pdf("Results/CorHeatmap.pdf", width=15, height = 15)
for(i in 1:length(meta_table)){
  i_col = meta_table[[i]]
  pheatmap(cor(ex), labels_col = i_col, labels_row = i_col, sub=names(meta_table)[i])
}
dev.off()


#### reduce dimension with PCA
ex.scale = t(scale(t(ex), scale = FALSE))
pc <- prcomp(ex.scale)
pdf("Results/PCA.pdf")
# by source name
pcr <- data.frame(pc$rotation[, 1:2], Group = meta_table[[4]])
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 3) + theme_bw()
# by Phenotype
pcr <- data.frame(pc$rotation[, 1:2], Group = meta_table[[5]])
ggplot(pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()


#### reduce dimension with MDS
ex.t <- t(ex)
d <- as.matrix(dist(ex.t))
fit <- cmdscale(d, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
pdf("Results/MDS.pdf")
# by source name
mds <- data.frame(x, y, Group = meta_table[[4]])
ggplot(mds, aes(x, y, color = Group)) + geom_point(size = 3) + theme_bw()
# by phenotype
mds <- data.frame(x, y, Group = meta_table[[5]])
ggplot(mds, aes(x, y, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()


#### reduce dimension with TSNE
tsne_results <- Rtsne(ex.t, perplexity = 10, check_duplicates = FALSE)
pdf("Results/TSNE.pdf")
# by source name
tsne <- data.frame(tsne_results$Y, Group = meta_table[[4]])
ggplot(tsne, aes(tsne_results$Y[, 1], tsne_results$Y[, 2], color = Group)) + geom_point(size = 3) + theme_bw()
# by phenotype
tsne <- data.frame(tsne_results$Y, Group = meta_table[[5]])
ggplot(tsne, aes(tsne_results$Y[, 1], tsne_results$Y[, 2], color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()


