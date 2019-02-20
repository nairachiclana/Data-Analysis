
options(warn=-1)

library(TCGAretriever)
library(RTCGAToolbox)
library(dplyr)
library(edgeR)
library(limma)

load(file="savedData/df_clinico.RData")

#load expression data 
getFirehoseDatasets()
brcaData <- getFirehoseData(dataset="BRCA", runDate="20160128",gistic2Date="20160128",forceDownload=F, clinical =TRUE, RNASeq2GeneNorm=TRUE)
#expression data
brca_rnaseq <- getData(brcaData,type = "RNASeq2GeneNorm")
#filter to keep only tumor data
brca_rnaseq.tumour <- brca_rnaseq[, which(as.numeric(substr(colnames(brca_rnaseq), 14,15)) < 10)]
#convert barcodes from the rnaseq dataset to sample barcodes, which identify the patients
colnames(brca_rnaseq.tumour) <- substr(colnames(brca_rnaseq.tumour), 1,12)
#remove duplicate samples
brca_rnaseq.tumour <- brca_rnaseq.tumour[, !duplicated(colnames(brca_rnaseq.tumour))]

#CANCER SUBTIPES

#Luminal A
luminal_samples <- df.tcga %>% dplyr::filter(`PAM50 mRNA` == "Luminal A")
luminal_barcodes <- luminal_samples$`Complete TCGA ID`
brca_rnaseq.luminal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% luminal_barcodes)]


#Triple Negativo
tnbc_samples <- df.tcga %>% dplyr::filter(`ER Status` == "Negative" & `PR Status` == "Negative" & `HER2 Final Status` == "Negative" & `PAM50 mRNA` != "Luminal A")
tnbc_barcodes <- tnbc_samples$`Complete TCGA ID`
brca_rnaseq.tnbc <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% tnbc_barcodes)]

df.rnaseq.all<-cbind(brca_rnaseq.luminal, brca_rnaseq.tnbc)
#numeric matrix of read counts
counts = df.rnaseq.all[apply(df.rnaseq.all,1,function(x) sum(x==0))<ncol(df.rnaseq.all)*0.8,]

#DESIGN MATRIX
df.luminal<-data_frame("sample"=colnames(brca_rnaseq.luminal), "status" = rep(0, length(colnames(brca_rnaseq.luminal))))
df.tnbc<-data_frame("sample"=colnames(brca_rnaseq.tnbc), "status" = rep(1, length(colnames(brca_rnaseq.tnbc))))
df.all<-rbind(df.luminal,df.tnbc)
design<-model.matrix(~ status, data=df.all)


#DEA MATRIX
dge<-DGEList(counts=counts)
A<-rowSums(dge$counts)
isexpr<- A > 100 
dge<-calcNormFactors(dge)

#transform RNA-Seq data for linear modelling (normalize)
v<-voom(dge[isexpr,], design, plot=FALSE)
#linear  model for each gene
fit <- lmFit(v, design)
#compute t-statistics by empirical Bayes moderation, generatin a precision weigth for each observation
fit <- eBayes(fit)
#extract top-ranked genes from the previous fit
diff.exp.df<-topTable(fit, coef = "status", n = Inf, sort = "p", p = 0.01) 
diff.exp.df$gene.name<-rownames(diff.exp.df)

save(diff.exp.df, counts, file="savedData/DEA.RData")

cat("Se ha realizado correctamente la matriz de expresión diferencial", "\n")


