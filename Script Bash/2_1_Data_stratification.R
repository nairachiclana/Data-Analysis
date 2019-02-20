options(warn=0)
load(file="savedData/clinical_expression.RData")

library(factoextra)

#AGE STRATIFICATION
colnames(clinical.expresion)[which(names(clinical.expresion) == "Age at Initial Pathologic Diagnosis")] <- "age"
clinical.expresion$age<-as.numeric(clinical.expresion$age)

data<-clinical.expresion
data$d1<-rep(1,length(clinical.expresion$age))
data.num=model.matrix(d1+age~age, data)

#clusters
data$cl= kmeans(data.num,2)$cluster
clinical.expresion$ClusteredAge<-cut(clinical.expresion$age, breaks=c(min(data$age[data$cl==1]),min(data$age[data$cl==2]),max(data$age[data$cl==2])))

cat("Age Stratification Done", "\n")

#CHANGE NAMES
colnames(clinical.expresion)[which(names(clinical.expresion) == "Metastasis-Coded")] <- "metastasis.coded"
colnames(clinical.expresion)[which(names(clinical.expresion) == "Node-Coded")] <- "node.coded"
colnames(clinical.expresion)[which(names(clinical.expresion) ==  "PAM50 mRNA")] <- "luminal.status"
colnames(clinical.expresion)[which(names(clinical.expresion) == "Tumor--T1 Coded")] <- "tumor.coded"
colnames(clinical.expresion)[which(names(clinical.expresion) == "Converted Stage")] <- "converged.stage"
colnames(clinical.expresion)[which(names(clinical.expresion) == "AJCC Stage")] <- "ajcc.stage"
colnames(clinical.expresion)[which(names(clinical.expresion) == "OS Time")] <- "time"
colnames(clinical.expresion)[which(names(clinical.expresion) == "OS event")] <- "status"
clinical.expresion$time<-as.numeric(clinical.expresion$time)
clinical.expresion$status<-as.numeric(clinical.expresion$status)
vbles.to.delete<-c("Complete TCGA ID","ER Status", "PR Status", "HER2 Final Status")
clinical.expresion<-clinical.expresion[, !(colnames(clinical.expresion) %in% vbles.to.delete)]

save(clinical.expresion, file="savedData/clinical_expression.RData")

