BiocManager::install("Biostrings")
BiocManager::install("limma")
BiocManager::install("illuminaHumanv3.db")
BiocManager::install("Biobase")
BiocManager::install("BiocGenerics")
BiocManager::install("BiocVersion")
BiocManager::install("Biostrings")
BiocManager::install("beadarray")
BiocManager::install("genefilter")
BiocManager::install("yarn")
BiocManager::install("arrayQualityMetrics")

install.packages("qusage")

install.packages("DESeq")
install.packages("gplots", dependencies = TRUE)
install.packages("dplyr")
install.packages("stringi")
install.packages("statmod")
install.packages("pheatmap")
install.packages("RColorBrewer", dependencies = TRUE)


library(pheatmap)
library(dplyr)
library(GEOquery)
library(beadarray)
library(affy)
library(limma)
library(genefilter)
library(yarn)
library(arrayQualityMetrics)
library(statmod)
library("DESeq")
library(gplots)
library(RColorBrewer)
library(qusage)

##Import the dataset using getGEO, accession number: GSE102451 - Illumina technology
mydata <- getGEO("GSE102451", GSEMatrix=TRUE, AnnotGPL=TRUE, getGPL=FALSE)
print(mydata)
View(mydata[[1]])

##The phenoData
View(pData(mydata[[1]])) 
##Matrix of the expression of each gene of each patient
View(exprs(mydata[[1]]))


arrayQualityMetrics(mydata[[1]], outdir= 'C:/Users/USER/Desktop/microArrays', force=TRUE)


##Extract the expression measurements.
##result: matrix with one row for each gene, and one column for each sample
head(exprs(mydata[[1]]))
ex.data <- exprs(mydata[[1]])
print(ex.data)
View(ex.data)    ##this is my expression matrix
summary(exprs(mydata[[1]]))   ##my data are not in logarithmic values
hist1 <- hist(ex.data) ##to see where my values are

##Examining the boxplot if the arrays are normalised.
bp1 <- boxplot(exprs(mydata[[1]]),outline=FALSE,axisnames=FALSE,las=2,col = c(rep("darkorange",27),"gold", "gold", rep("darkolivegreen1",29), "cornsilk4") )  #normalized dataset

##M-A plots:compare the red and green channels from a two-colour microarray
##log.it = TRUE uses logarithm of the values
mva.pairs(ex.data[,1:5],log.it = TRUE)

bp1_log <- boxplot(log2(ex.data),outline=FALSE,las=2)
hist1 <- hist(log2(ex.data))
##Normalise using normaliseIllumina     
##mydata.norm <- normaliseIllumina(mydata[[1]])  ##NOT DO IT, THE DATA IS ALREADY NORMALISED??
##View(exprs(mydata.norm))
mydata.norm.log <- log2(ex.data)
View(mydata.norm.log)

bp2_log <- boxplot(mydata.norm.log,outline=FALSE)

mva.pairs(mydata.norm.log[,1:4])

limma::plotMA(ex.data) 
limma::plotMA(mydata.norm.log)

##The annotation for the features
##all the rows of this feature matrix are in the same order as the expression matrix
#fData(mydata[[1]])[1:5,1:5]           ##ERROR because my fData is empty
##The fdata=phenoData has info on each gene (metadata)          
#all(rownames(fData(mydata)) == rownames(exprs(mydata)))  

##Dataset information from >featureData we prnt the >data
pData(mydata[[1]])[1:5,]
colnames(pData(mydata[[1]]))
View(pData(mydata[[1]]))
##I have 3 sample groups:
###HCC tissue from cirrhotic liver [1:27]
###Pre-malignant dysplastic nodule from cirrhotic liver [28:29]
###Non-tumor cirrhotic liver [30:59]

##Differential expression analysis

##In order to do paired t test, remove the last patient that does not have a pair
subset <- mydata[[1]][,-c(59)]
ex.subset <- exprs(subset)
dim(ex.subset)  ##20761    58
ex.subset <- unique(ex.subset)
dim(ex.subset)  ##20761    58


View(subset)

pd <- pData(subset)
View(pd)
SampleGroup <- pd$source_name_ch1
View(SampleGroup)
subset.norm.log <- log2(ex.subset)
View(subset.norm.log)

##Filtering before statistic analysis
##filter the data so only the top 50% most-variable genes get analysed
##increase our power to detect differential expression
##mydata.expFilt <- varFilter(subset)  
##View(mydata.expFilt)      ##Because of the Bayes step do not filter

##HCC and pre-malignant -> 0, Non-tumor cirrhotic -> 1
trt <- factor(rep(0:1, each = 29))

pairs <- factor(rep(1:29, 2))

##Correlation based on paires -> "Tumor","Cirrhosis"
design <- model.matrix(~0+trt)

colnames(design) <- c("Tumor","Cirrhosis")

##Correlations given my design on the "Tumor","Cirrhosis" paires
corfit <- duplicateCorrelation(subset.norm.log,design,block=pairs)
corfit$consensus

##The lmFit funcion is used to fit the model to the data. 
fit <- lmFit(subset.norm.log, design,block=pairs,correlation=corfit$consensus)

contrasts <- makeContrasts(Tumor-Cirrhosis,levels=design)

fitc <- contrasts.fit(fit,contrasts)

fit2 <- eBayes(fitc)
##Paired t-test for my expression data
##Sort by the absolute value of logFC
topTable(fit2)
volcanoplot(fit2)

volcanoplot(fit2,highlight=5,names = row.names(fit2))
with(subset(genes.table, adj.P.Val<.05 & abs(logFC)>1), points(logFC, -log10(P.Value), pch=20, col="cyan2"))

##All genes sorted by absolute value of logFC
genes.table <- topTable(fit2, number = length(ex.data))

#Topgenes - filter my output (p-value < 0.05 and (logFC > 1 or logFC < -1))
##20761

topgenes <- subset(genes.table,genes.table$P.Value < 0.05 & (genes.table$logFC > 1 | genes.table$logFC < -1))

##208

mytopGenes <- row.names(topgenes)
View(mytopGenes)

write.table(mytopGenes, quote = FALSE, "C:/Users/USER/Documents/Biology/microarrays/mytopGenes.txt", sep="\t", row.names = FALSE, col.names = FALSE)

##Get the expressions of my top genes
###Gives a subset dataframe with my top genes
exp.mytopgenes <- subset.norm.log[mytopGenes, ] ##Log or not log data
View(exp.mytopgenes)
     
write.table(exp.mytopgenes, file="C:/Users/USER/Documents/Biology/microarrays/exp_topGenes.csv", sep=",", quote=F, row.names=T)

##Needs DESeq and featmap libraries
pheatmap(exp.mytopgenes)

#heatmap_xaxis <- c('P_1','P_2','P_3','P_4','P_5','P_6','P_7','P_8','P_9','P_10','P_11','P_12','P_13','P_14','P_15','P_16','P_17','P_18','P_19','P_20','P_21','P_22','P_23','P_24','P_25','P_26','P_27','P_28','P_29','P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12','P13','P14','P15','P16','P17','P18','P19','P20','P21','P22','P23','P24','P25','P26','P27','P28','P29')

labels <- c(rep("Tumor",29),rep("Cirrhosis",29))

my_sample_col <- data.frame(labels)
row.names(my_sample_col) <- colnames(exp.mytopgenes)

pheatmap(exp.mytopgenes, annotation_col = my_sample_col, show_rownames = FALSE,fontsize = 10, fontsize_col = 5.4)  


# #HeatMap
# rmeans <- rowMeans(subset.norm.log)
# subset.norm.log.sort <- subset.norm.log[order(rmeans,decreasing=T),]
# heatmap(subset.norm.log.sort[1:500,])



##Perform hierarchical clustering to obtain gene clusters. 
##I use the excellent dendextend to plot a simple dendrogram.

my_hclust_gene <- hclust(dist(exp.mytopgenes), method = "complete")

# install
install.packages("dendextend")
library(dendextend)

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

###Enrichment analysis
### A pairVector containing information on which samples should be paired together.
labels <- c(rep("Tumor",29),rep("Cirrhosis",29))
print(labels)

##single character string: how you want to compare the groups in your data
contrast <- "Tumor - Cirrhosis"

##A vector describing a single gene set
gene.names <- c(row.names(subset.norm.log))
gene.names

##The qusage method returns a QSarray object containing statistical data on both the
##distributions of individual genes and on the pathway itself
##For my top genes only
qs.results <- qusage(exp.mytopgenes, labels, contrast, gene.names)

##To find out if our geneSet is significantly enriched for this comparison:
pdf.pVal(qs.results)
##0.0001093849
plot(qs.results, xlab="Gene Set Activation")

pairs <- c(1:29,1:29)


MSIG.geneSets <- read.gmt("C:/Users/USER/Documents/Biology/microarrays/c2.cp.kegg.v7.0.symbols.gmt")

summary(MSIG.geneSets[1:5])
MSIG.geneSets[2]
qs.results.msig <- qusage(exp.mytopgenes, labels, contrast, MSIG.geneSets,pairVector=pairs)

numPathways(qs.results.msig)
p.vals = pdf.pVal(qs.results.msig)
head(p.vals)
q.vals = p.adjust(p.vals, method="fdr")
head(q.vals)
plot(qs.results.msig,cex=0.6)

#######################################
##Diagnostic checks for top genes
dotchart(subset["VIPR1",])
boxplot(subset["VIPR1",]~SampleGroup,las=2)