options(warn=-1)

library(pROC)
library(caret)
library(glmnet)

load(file="savedData/replicated_df.RData")

args = commandArgs(trailingOnly=TRUE)
k = as.numeric(args[1])

cat("La k recibida es", k)

total=10*k
pb <- txtProgressBar(min = 0, max=total, style = 3)

#LASSO AND RIDGE
cf<-createFolds(replicated.df$status, k=k)
l.auc<-vector()
l.alpha<-vector()
l.lambda<-vector()

for(j in 1:k) {
  test <-replicated.df[unlist(cf[j]),]
  train <-replicated.df[-unlist(cf[j]),]
  
  train.matrix<-model.matrix(status~.,train)
  test.matrix<-model.matrix(status~.,test)
  status.train<-train$status
  status.test <- test$status
  
  lista<-list()
  aucs<-c()
  lambdas<-c()
  
  init.lr10<-Sys.time()
  
  for (i in 0:10) { 
    
    lambdas<-c(lambdas,cv.glmnet(train.matrix, status.train, family = "binomial",
                                 nfold = 10, type.measure = "auc", paralle = TRUE, alpha = 1)$lambda.1se)
    lista[[(i+1)]]<-glmnet(train.matrix, status.train, family = "binomial",lambda = lambdas[i+1],alpha = i/10)
    
    aucs<-c(aucs,roc(status.test, as.numeric(predict(lista[[(i+1)]], test.matrix, type = "response"))))
    
    setTxtProgressBar(pb, i*j)
  }
  
  #train 
  lasso.model<-glmnet(train.matrix, status.train, family = "binomial",lambda = min(lambdas),alpha = (which.min(lambdas)-1)/10) 
  #prediction
  lasso.prob.1 <- predict(lasso.model,type="response",newx = test.matrix) 
  new.auc<-roc(status.test,as.numeric(lasso.prob.1))$auc 
  
  l.auc<-c(l.auc, new.auc)
  l.alpha<-c(l.alpha, (which.min(lambdas)-1)/10)
  l.lambda<-c(l.lambda, min(lambdas))
  
}

end.lr10<-Sys.time()
time.lr10=end.lr10-init.lr10




cat("\n", "Los alphas de las iteraciones son: ", l.alpha, "\n")
cat("Los lambda de las iteraciones son: ", l.lambda, "\n")
cat("Los auc de las iteraciones son: ", l.auc, "\n")
cat("El tiempo de ejecución ha sido: ", time.lr10, "\n")


