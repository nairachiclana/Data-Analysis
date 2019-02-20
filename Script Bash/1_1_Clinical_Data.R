options(warn=-1)

library(readxl)

sample_data<-read_excel("nature11412-s2/Supplementary Tables 1-4.xls", sheet = 1, skip = 1)
df.tcga<-as.data.frame(sample_data)

#vbles to keep
variables_to_keep<-c("Complete TCGA ID","Gender" , "Age at Initial Pathologic Diagnosis", "ER Status", "PR Status",  "HER2 Final Status","Tumor","Tumor--T1 Coded" , "Node", "Node-Coded", "Metastasis", "Metastasis-Coded", "AJCC Stage", "Converted Stage","OS event" ,"OS Time" , "PAM50 mRNA")

#eliminate rows from colum list
eliminate.na.rows<-function(df,variables_need_imputation) {
  for(variable in variables_need_imputation) {
    col<-df[,variable]
    df<-df[col!="NA",]
  }
  return(df)
}

df.tcga<-df.tcga[ , (names(df.tcga) %in% variables_to_keep)]
df.tcga.imp<-eliminate.na.rows(df.tcga, variables_to_keep)
df.tcga$"OS Time"<-as.numeric(df.tcga$"OS Time")
df.tcga$"OS event"<-as.numeric(df.tcga$"OS event")



cat("Dimensiones del dataset tras la imputación: ", dim(df.tcga.imp), "\n")

df.tcga<-df.tcga.imp

save(df.tcga, file="savedData/df_clinico.RData")