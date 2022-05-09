install.packages("BiocManager")
install.packages("TCGAbiolinks")
install.packages("DESeq2")
BiocManager::install("PCAtools")
BiocManager::install('biomaRt')
BiocManager::install("apeglm")
BiocManager::install("edgeR")
library("BiocManager")
library('TCGAbiolinks')
library("DESeq2") 
library("biomaRt") 
library("apeglm") 
library("PCAtools")
library('edgeR')

TCGAdata <- GDCquery(
  project = 'TCGA-BLCA',
  data.category = 'Transcriptome Profiling',
  experimental.strategy = 'RNA-Seq',
  workflow.type = 'HTSeq - Counts')

GDCdownload(query = TCGAdata)

bcdata <- GDCprepare(query = TCGAdata, save = TRUE, save.filename = 'expresion.rda')


rna <- as.data.frame(SummarizedExperiment::assay(data))
clinical <- data.frame(data@colData)
table(clinical$definition)
table(substr(colnames(rna),14,14))

clinical$paper_mRNA.cluster <-  gsub(" ", "_", clinical$paper_mRNA.cluster)

clinical$paper_mRNA.cluster <- as.factor(clinical$paper_mRNA.cluster)

levels(clinical$paper_mRNA.cluster)
clinical$paper_mRNA.cluster <- relevel(clinical$paper_mRNA.cluster, ref = ".")
clini_new = clinical[!is.na(clinical$paper_mRNA.cluster), ]
rna_new = rna[ ,!is.na(clinical$paper_mRNA.cluster)]

dds <- DESeqDataSetFromMatrix(countData = rna_new,
                              colData = clini_new,
                              design = ~ paper_mRNA.cluster)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


dds <- DESeq(dds) 

res <- results(dds, alpha = 0.05,  altHypothesis = "greaterAbs", lfcThreshold = 1.5) # alpha controls FDR rate


resLFC <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
resLFC.Ordered<-resLFC[with(resLFC, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]


ens2symbol<-function(ids){
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- getBM(filters= "ensembl_gene_id", 
                 attributes= c("ensembl_gene_id","hgnc_symbol"),
                 values=ids, mart= mart)
  return(genes)
}

df <- ens2symbol(row.names(res))

res_df <- as.data.frame(res)                 
res_df$ensembl_gene_id <- row.names(res_df)
res_df <- merge(df,res_df, by = "ensembl_gene_id")
resOrdered<-res_df[with(res_df, order(abs(log2FoldChange), padj, decreasing = TRUE)), ]

write.csv(resOrdered, "KMT2D_dfnew.csv")

ylim <- c(-6.5,6.5)
drawLines <- function() abline(h=c(-1.5,1.5),col="dodgerblue",lwd=2)
plotMA(resLFC, ylim=ylim); drawLines()
