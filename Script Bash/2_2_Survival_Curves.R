options(warn=0)
load(file="savedData/clinical_expression.RData")

library(KMsurv)
library(survMisc)
library(survminer)
library(flexsurv)
library(ggfortify)

#Global
gg.global<-ggsurvplot(survfit(Surv(time,status)~1, data=clinical.expresion, type="kaplan-meier"),  xlab = "Tiempo", censor = T, conf.int=T, ylab = "Survival Probability", title = "Survival probability",  ggtheme = theme_bw())

#Clustered Age
gg.age<-survfit(Surv(time,status)~ClusteredAge, clinical.expresion, conf.type = "log-log") %>% ggsurvplot(title = "Supervivencia por grupo edad",  legend.title = "Grupo edad", pval=TRUE)

#Metastasis.coded
gg.metas.cod<-survfit(Surv(time,status)~metastasis.coded, clinical.expresion, conf.type = "log-log") %>% ggsurvplot(title = "Supervivencia por estado metastasis",  legend.title = "Estado metastasis", pval=TRUE,  ggtheme = theme_bw())

#Node.coded
gg.node.cod<-survfit(Surv(time,status)~node.coded, clinical.expresion, conf.type = "log-log") %>% ggsurvplot(title = "Supervivencia por Codificación de tumor",  legend.title = "Tumor", pval=TRUE,  ggtheme = theme_bw())

#Node
gg.metas.cod<-survfit(Surv(time,status)~Node, clinical.expresion, conf.type = "log-log") %>% ggsurvplot(title = "Supervivencia por Nodos",  legend.title = "Tumor", pval=TRUE,  ggtheme = theme_bw())

#Tumor
gg.tumor<-survfit(Surv(time,status)~Tumor, clinical.expresion, conf.type = "log-log") %>% ggsurvplot(title = "Supervivencia por Tumor",  legend.title = "Tumor", pval=TRUE,  ggtheme = theme_bw())

#Luminal status
gg.luminal<-survfit(Surv(time,status)~luminal.status, clinical.expresion, conf.type = "log-log") %>% ggsurvplot(title = "Supervivencia por subtipo", legend.title = "Subtipo", pval=TRUE,  ggtheme = theme_bw())



gg.list<-list()
  
gg.list[[1]]<-gg.global
gg.list[[2]]<-gg.age
gg.list[[3]]<-gg.metas.cod
gg.list[[4]]<-gg.node.cod
gg.list[[5]]<-gg.tumor
gg.list[[6]]<-gg.luminal


ggsave("results/survival_curves.png", arrange_ggsurvplots(gg.list, ncol=2, nrow=3),device = 'png', width = 13, height = 10)

