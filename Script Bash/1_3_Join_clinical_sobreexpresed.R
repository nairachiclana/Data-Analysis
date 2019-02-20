options(warn=-1)

library(dplyr)

load(file="savedData/DEA.RData")
load(file="savedData/df_clinico.RData")


args = commandArgs(trailingOnly=TRUE)
lfc = as.numeric(args[1])
pvalue = as.numeric(args[2])

cat("El pvalue introducido es", pvalue, "el lfc introducido es: ", lfc, "\n")
#MATRIX SOBREEXPRESADOS

#filter by lfc
sobreexpresed.genes<-rbind(diff.exp.df[diff.exp.df$logFC<(-lfc),] , diff.exp.df[diff.exp.df$logFC>lfc,])
#filter by pvalue
sobreexpresed.genes<-sobreexpresed.genes[ sobreexpresed.genes$adj.P.Val<pvalue,]

clinical.data<-df.tcga

#JOIN
transpuesta<-data.frame(t(counts[sobreexpresed.genes$gene.name,]))
transpuesta$`Complete TCGA ID` <- rownames(transpuesta)

compact<- clinical.data[c(clinical.data$`Complete TCGA ID`) %in% transpuesta$`Complete TCGA ID`,]
clinical.expresion<-inner_join(compact, transpuesta)

save(clinical.expresion, file="savedData/clinical_expression.RData")

cat("Las dimensiones de la matriz unida son: ", dim(clinical.expresion), "\n")
