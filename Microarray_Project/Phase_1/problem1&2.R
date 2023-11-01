library(GEOquery)
library(limma)
library(Biobase)
library(pheatmap)
library(reshape2)
library(plyr)
library(ggplot2)
library(stringr)

##########
seriesName <- "GSE48558"
platformName <- "GPL6244"

gset <- getGEO(seriesName,GSEMatrix = TRUE, AnnotGPL = TRUE , destdir = "D:/Term/term7/BIO/project")

if (length(gset) > 1){
  idx <- grep(platformName, attr(gset, "names"))
} else {
  idx <- 1
}
gset <- gset[[idx]]

gset<- gset[,which(gset$source_name_ch1 == "AML Patient" | gset$`phenotype:ch1` == "Normal")]


func <- function(x) {
  if (gset$source_name_ch1[x] == "AML Patient") {
    return("Test")
  } else {
    spll <- strsplit2(gset$source_name_ch1[x] , "\+")[1, 1]
    return(paste0("Normal_" , spll))
  }
}



gr <- sapply(1:length(gset$`phenotype:ch1`) , func)

###########

expr <- exprs(gset)
print(paste0(max(expr)))
print(paste0(min(expr)))

##########
expr <- log2(1 + expr)
exprs(gset)<- expr
pdf("D:/Term/term7/BIO/project/plot/boxplot.pdf" , width = 32)
boxplot(expr)
dev.off()
#######
expr <- normalizeQuantiles(expr)
exprs(gset) <- expr
pdf("D:/Term/term7/BIO/project/plot/boxplot_afterNormalize.pdf" , width = 32)
boxplot(expr)
dev.off()