setwd("C:/Users/Mohammad/Desktop/Project/phase2/R-wd/")

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
ex <- exprs(gset)


#### drawing box plot for quality control
pdf("Results/boxplot.pdf",width = 64)
boxplot(ex)
dev.off()


#### reduce dimension with PCA
ex.scale = t(scale(t(ex), scale = FALSE))
pc <- prcomp(ex.scale)


### filtering ignored data
pcr <- data.frame(pc$rotation[, 1:2], Group = meta_table[[4]])
included_groups = c("CD34+HSPC", "Monocytes", "AML Patient")
included_meta = meta_table[meta_table$Source.name %in% included_groups,]
ex.included <- ex[, which(meta_table[[4]] %in% included_groups != 0)]

### draw pca
pdf("Results/PCA.pdf")
# by source name
included_pcr = pcr[pcr$Group %in% included_groups,]
ggplot(included_pcr, aes(PC1, PC2, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()


#### reduce dimension with MDS
ex.t <- t(ex.included)
d <- as.matrix(dist(ex.t))
fit <- cmdscale(d, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
pdf("Results/MDS.pdf")
# by source name
mds <- data.frame(x, y, Group = included_meta[[4]])
ggplot(mds, aes(x, y, color = Group)) + geom_point(size = 3) + theme_bw()
# by phenotype
mds <- data.frame(x, y, Group = included_meta[[5]])
ggplot(mds, aes(x, y, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()


#### reduce dimension with TSNE
tsne_results <- Rtsne(ex.t, perplexity = 8, check_duplicates = FALSE)
pdf("Results/TSNE.pdf")
# by source name
tsne <- data.frame(tsne_results$Y, Group = included_meta[[4]])
ggplot(tsne, aes(tsne_results$Y[, 1], tsne_results$Y[, 2], color = Group)) + geom_point(size = 3) + theme_bw()
# by phenotype
tsne <- data.frame(tsne_results$Y, Group = included_meta[[5]])
ggplot(tsne, aes(tsne_results$Y[, 1], tsne_results$Y[, 2], color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()



### making table
contrast_source_names = meta_table[[4]]
for(i in 1:length(source_names)){
  contrast_source_names[i] <- gsub(" ", "", contrast_source_names[i])
  contrast_source_names[i] <- gsub("Monocytes", "Health", contrast_source_names[i])
  contrast_source_names[i] <- gsub("CD34\\+HSPC", "Health", contrast_source_names[i])
}
gset$group <- contrast_source_names
design <- model.matrix(~group + 0, gset)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(contrasts='groupAMLPatient-groupHealth', levels=design)
contrast_fit <- contrasts.fit(fit, cont.matrix)
contrast_fit <- eBayes(contrast_fit, 0.01)
top_table <- topTable(contrast_fit, adjust="fdr", sort.by="B", number=Inf)
top_table <- subset(top_table, select=c("Gene.symbol", "Gene.ID", "adj.P.Val","logFC"))
head(top_table)


### get top of table and high frequented in Aml
aml.up <- subset(top_table, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file="Results/highAml.txt", row.names=FALSE, col.names = FALSE, quote = FALSE, sep="\t")

### get tail of table and low frequented in Aml
aml.low <- subset(tT, logFC <- 1 & adj.P.Val < 0.05)
aml.low.genes <- unique(as.character(strsplit2(aml.low$Gene.symbol, "///")))
write.table(aml.low.genes, file="Results/lowAml.txt", row.names=FALSE, col.names = FALSE, quote = FALSE, sep="\t")
