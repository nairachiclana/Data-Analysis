options(warn=-1)
load(file="savedData/clinical_expression.RData")

library(FSelector)
library(survival)

args = commandArgs(trailingOnly=TRUE)
filter.subset.size = as.numeric(args[1])
cox.subset.size = as.numeric(args[2])

#CLINICAL:
#Estas variables han sido elejidas a mano sin ningún procedimiento, solamente por observación de los datos
clinical.expresion<-clinical.expresion[ , !(names(clinical.expresion) %in% "tumor.coded")]
clinical.expresion<-clinical.expresion[ , !(names(clinical.expresion) %in% "Metastasis")]
clinical.expresion<-clinical.expresion[ , !(names(clinical.expresion) %in% "age")]
clinicas.seleccionadas<-names(clinical.expresion)[1:10]
cat("Las variables clinicas finales son: ",clinicas.seleccionadas, "\n")

#EXHAUSTIVE SEARCH
all<-names(clinical.expresion)
'%ni%'<-Negate('%in%')
#with genes
vbles<-all[(all %ni% c("time", "status"))]


pvalue.acumulado<-list()

for(nivel in 1:length(vbles)) { #cada nivel (#vbles=#nivel)
  best.pvalue<-1
  for(v in vbles) { #vbles de cada nivel
    ifelse(nivel==1, att<-v, att<-c(v, lista.acumulada))
    ec <-as.simple.formula(att,"Surv(time,status)")
    sum<-summary(coxph(ec, data=clinical.expresion))
    new.p<-sum$logtest["pvalue"]
    if(new.p <0.05 && new.p<best.pvalue) {
      best.pvalue<-new.p
      best.vble.nivel<-v
    }
  }
  
  vbles<-vbles[vbles != best.vble.nivel]
  
  ifelse(nivel!=1,  pvalue.acumulado<-c(pvalue.acumulado, best.pvalue), first.pvalue<-best.pvalue)
  ifelse(nivel==1, lista.acumulada<-best.vble.nivel,lista.acumulada<-c(lista.acumulada, best.vble.nivel))
  ifelse(nivel==1, p.acum.ac<-best.pvalue, p.acum.ac<-pvalue.acumulado[nivel-1])
  if(length(lista.acumulada)!=nivel || p.acum.ac>best.pvalue[1] || length(lista.acumulada)==cox.subset.size) break
}

cat("Variables resultantes de Busqueda exhasustiva con Cox:",  lista.acumulada, "\n")

#FILTERS
formula<-as.formula(paste("status",".",sep="~"))
#Chi-squared
chi.weigths<-chi.squared(formula, clinical.expresion)
best.chi.weigths<-cutoff.k(chi.weigths,filter.subset.size)
cat("Las variables filtradas por chi-squared son: ",best.chi.weigths, "\n")
#Forest-importance
forest.weigths<-random.forest.importance(formula, clinical.expresion)
best.forest.weigths<-cutoff.k(forest.weigths,filter.subset.size)
cat("Las variables filtradas por random forest son: ",best.forest.weigths, "\n")

coinc1<-best.chi.weigths[best.chi.weigths %in% forest.weigths]
coinc2<-best.forest.weigths[best.forest.weigths %in% best.chi.weigths]
coincidentes<-unique(c(coinc1,coinc2))
coincidentes.con.cox<-lista.acumulada[lista.acumulada %in% coincidentes]

cat("Las variables coincidentes entre filters son: ", coincidentes, "\n")
cat("Las variables coincidentes entre filters y cox son: ", coincidentes.con.cox, "\n")

#Hemos visto un ejemplo de como funciona, pero por seguridad de que los parametros de tamaño seleccinados no sean correctos para este caso, 
#guardaremos a mano el dataset final procesado y explicado en el archivo html

final.vbles<-c("Gender", "Tumor", "Node", "node.coded", "metastasis.coded", "ajcc.stage", "converged.stage", "status", "time", "luminal.status", "SPDEF", "PRR15", "TFF3", "PGLYRP2","ERBB4", "SYTL5","AR","TFF1", "CST9", "KLHDC7A", "CEACAM5", "A2ML1","GABBR2", "GABRP", "KRT16","KLK6")
filtered.df<-clinical.expresion[ , (names(clinical.expresion)) %in%  final.vbles]
cat("VARIABLES FINALES:", names(filtered.df), "\n")

save(clinical.expresion, file="savedData/clinical_expression.RData")
save(filtered.df,file="savedData/filtered_df.RData" )