options(warn=-1)


args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])

load(file="savedData/replicated_df.RData")

library(caret)
library(stats) #glm
library(nnet) #nnet
library(rpart) #dt
library(grDevices)
library(pROC)

#AUC FUNCTION
make.auc <-function(mod.estimado, datos) {
  if(is.factor(mod.estimado)) mod.estimado<-ifelse(mod.estimado=="X1", 1,0)
  roc.curve <-roc(datos$status, mod.estimado, smooth=FALSE, auc=TRUE)
  auc.result <-roc.curve$auc
  return(auc.result)
}

compute.models <-function(train.data, test.data, model, formula) {
  
  #Decision Trees
  if(model=="dt") {
    modelo<-rpart(formula,data=train.data,control=rpart.control(minsplit=1, cp=0.0059));
    mod.estimado<-predict(modelo, newdata=test.data, type="matrix");
    mod.estimado<-as.matrix(mod.estimado)
  }
  
  #Regresión logística
  if(model=="glm") {
    modelo<-glm(formula, train.data, family=binomial("logit"));
    mod.estimado<-predict(modelo, newdata=test.data, type="response");
  }
  
  #Neural network
  else if(model=="nnet") { 
    modelo<-nnet(formula,data=train.data, size=5, decay=0, maxit=7, trace=F);
    mod.estimado<-predict(modelo, newdata=test.data, type="raw");
  }
  
  if(is.factor(mod.estimado)) mod.estimado<-ifelse(mod.estimado=="X1", 1, 0)
  ifelse(model=="dt", roc.test<-roc(as.numeric(test.data$status), as.numeric(mod.estimado[,dim(mod.estimado)[2]]), smooth=F, auc=T), roc.test<-roc(as.numeric(test.data$status), mod.estimado, smooth=F, auc=T))
  auc.test <-roc.test$auc
  return(list(auc=auc.test,roc=roc.test)) 
}

load(file="savedData/replicated_df.RData")

pb <- txtProgressBar(min = 0, max = k, style = 3)

  auc.glm<-list()
  auc.dt<-list()
  auc.nnet<-list()
  cf<-createFolds(replicated.df$status, k=k)
  formula=as.formula(paste("status", ".", sep="~"))
  
  init.cv<-Sys.time()
  
  for(i in 1:k) {
    
    test <-replicated.df[unlist(cf[i]),]
    train <-replicated.df[-unlist(cf[i]),]
    train$status <-as.factor(train$status)
    test$status <-as.factor(test$status)
    
    auc.glm<-c(auc.glm,compute.models(train, test, "glm", formula)$auc);
    auc.dt<-c(auc.dt,compute.models(train, test, "dt", formula)$auc);
    auc.nnet<-c(auc.nnet,compute.models(train, test, "nnet", formula)$auc);
    
    setTxtProgressBar(pb, k)
  }
  
  end.cv<-Sys.time()
  
  #means
  mean.glm=mean(unlist(auc.glm))
  mean.dt=mean(unlist(auc.dt))
  mean.nnet=mean(unlist(auc.nnet))
  means.cv<-list(mean.glm, mean.dt, mean.nnet)
  #boxplot   
  png("results/predictors_auc_boxplot.png")
  boxplot(as.numeric(auc.glm),as.numeric(auc.dt),as.numeric(auc.nnet), names=c("GLM", "DT", "NNET"))
  names=c("GLM", "DT", "NNET")
  points(x=1:3,y=means.cv,pch=19,col="green4")
  text(x=1:3,y=means.cv, format(round(as.numeric(means.cv), 2), nsmall = 3), col="mediumseagreen", pos=1)
  dev.off()
  
cat("El auc de glm ha sido: ", mean.glm, "\n")
cat("El auc de dt ha sido: ", mean.dt, "\n")
cat("El auc de nnet ha sido: ", mean.nnet, "\n")
cat("El tiempo de ejecución ha sido: ", end.cv-init.cv, "\n")


