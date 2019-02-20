---
title: 'P3: Survival Analysis'
author: "Naira María Chiclana García 44717497T"
date: "January 2019"
output: 
  html_document:
    self_contained: false
    keep_md: true
    toc: true
    toc_float: true
---


Objective: **Use the Cox regression model and ANN for making survival predictions**


 





```r
library(survival)
library(survminer)
library(survMisc)
library(ggplot2)
```



# Data set

Lung cancer data from NCCTG (North Central Cancer Treatment Group)


```r
df.lung<-lung
head(df.lung)
```

```
##   inst time status age sex ph.ecog ph.karno pat.karno meal.cal wt.loss
## 1    3  306      2  74   1       1       90       100     1175      NA
## 2    3  455      2  68   1       0       90        90     1225      15
## 3    3 1010      1  56   1       0       90        90       NA      15
## 4    5  210      2  57   1       1       90        60     1150      11
## 5    1  883      2  60   1       0      100        90       NA       0
## 6   12 1022      1  74   1       1       50        80      513       0
```

Colums:

- `inst:` Intitution code.

- `time`:	Survival time in days.

- `status`:	censoring status 1=censored, 2=dead. The censored data are those patients in which the event haven't occur by the time of data collection, this data have also to be taken into account.

- `age`:	Age in years.

- `sex`:	Male=1 Female=2.

- `ph.ecog`:	ECOG performance score (0=good 5=dead).

- `ph.karno`:	Karnofsky performance score (bad=0-good=100) rated by physician.

- `pat.karno`:	Karnofsky performance score as rated by patient.

- `meal.cal`:	Calories consumed at meals.

- `wt.loss`:	Weight loss in last six months.

**ECOG performance: ** Scale developed by Eastern Cooperative Oncology Group describes the patient's level of functioning in terms of theirr abiity to care for themself, daily activity and physical ability. Being: 0 fully active, 1 Restricted in physically strenuous activity but ambulatory and able to carry out work of a light or sedentary nature, 2 Ambulatory and capable of all selfcare but unable to carry out any work activities, 3 Capable of only limited selfcare, 4 Completely disabled and 5 dead.

**Karnofsky performance: ** As ECOG performance, allows patients to be classified as to their functional impairment. It's somehow "inverse" to ECOG, Being Karno=100/90 equivalent to ECOG=0 (normal and able to carry on normal activity) and Karno=0 equivalent to ECOG=5 (Dead).



### Stratification of  covariates

#### Age
 
We're going to use clustering to divide the age in groups.


```r
library(factoextra)
```

```
## Welcome! Related Books: `Practical Guide To Cluster Analysis in R` at https://goo.gl/13EFCZ
```

```r
data<-df.lung
data$d1<-rep(1,length(df.lung$age))
#Convert factor variables to numeric 'dummy' variables
data.num=model.matrix(d1+age~age, data)
#número óptimo de clusters con su puntuación silohuette
fviz_nbclust(data.frame(data$age,data$d1), kmeans, method = "silhouette") +geom_vline(xintercept = 2, linetype = 2)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

The optimal number of clusters for the age covariate is 2, with an average silohuette index of 0.6 *(Saw in Practice 1)*. We're better going to make three clusters since we want more differenciate groups and the silhouette index for this k is also good.


```r
#Assign clusters
data$cl= kmeans(data.num,3)$cluster
#Distribution of data
table(data$cl)
```

```
## 
##  1  2  3 
## 85 94 49
```

```r
df.lung$ClusteredAge<-cut(df.lung$age, breaks=c(min(data$age[data$cl==2]),min(data$age[data$cl==3]),min(data$age[data$cl==1]),max(data$age[data$cl==1])))
table(df.lung$ClusteredAge)
```

```
## 
## (39,56] (56,67] (67,82] 
##      56      93      77
```

The third cluster have too little number of clusters, we're going to keep only two.


```r
data$cl= kmeans(data.num,2)$cluster
table(data$cl)
```

```
## 
##   1   2 
## 129  99
```

```r
table(data$age[data$cl==1])
```

```
## 
## 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 80 81 82 
##  7 11 11  8  7  8 10 11 10  7  7  6 10  5  5  2  2  1  1
```

```r
table(data$age[data$cl==2])
```

```
## 
## 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 
##  2  1  1  1  1  5  1  1  1  4  2  6  2  2  9  4  6  9  9  8  8 11  5
```

```r
df.lung$ClusteredAge<-cut(df.lung$age, breaks=c(min(data$age[data$cl==2]),min(data$age[data$cl==1]),max(data$age[data$cl==1])))
table(df.lung$ClusteredAge)
```

```
## 
## (39,62] (62,82] 
##     104     122
```
The distribution seems to be more balanced now.


#### Karnofsky scores

According to frequently medical separation groups:

- $[100, 80]$: Able to carry on normal activity and to work; no special care needed.

- $[70,50]$: Unable to work; able to live at home and care for most personal needs; varying amount of assistance needed.

- $[40, 0]$: Unable to care for self; requires equivalent of institutional or hospital care; disease may be progressing rapidly.


```r
df.lung$group.ph.karno<-cut(df.lung$ph.karno, breaks=c(0,40,70,100))
table(df.lung$group.ph.karno)
```

```
## 
##   (0,40]  (40,70] (70,100] 
##        0       57      170
```

There's no patients unables to taking care of themselves accordingo to physician.


```r
df.lung$group.pat.karno<-cut(df.lung$pat.karno, breaks=c(0,40,70,100))
table(df.lung$group.pat.karno)
```

```
## 
##   (0,40]  (40,70] (70,100] 
##        4       75      146
```

According to patients, there are some of them.

The dataframe with the new columns added:



```r
head(df.lung)
```

```
##   inst time status age sex ph.ecog ph.karno pat.karno meal.cal wt.loss
## 1    3  306      2  74   1       1       90       100     1175      NA
## 2    3  455      2  68   1       0       90        90     1225      15
## 3    3 1010      1  56   1       0       90        90       NA      15
## 4    5  210      2  57   1       1       90        60     1150      11
## 5    1  883      2  60   1       0      100        90       NA       0
## 6   12 1022      1  74   1       1       50        80      513       0
##   ClusteredAge group.ph.karno group.pat.karno
## 1      (62,82]       (70,100]        (70,100]
## 2      (62,82]       (70,100]        (70,100]
## 3      (39,62]       (70,100]        (70,100]
## 4      (39,62]       (70,100]         (40,70]
## 5      (39,62]       (70,100]        (70,100]
## 6      (62,82]        (40,70]        (70,100]
```

### Data imputation 

As  we will see in [Selection of best covariates for cox model](#linkcovariatesselection), the used variables will be *Clustered age* and the groups of *Karnofsky scores*. We will imputate the rows in  which any of these values is missing.


```r
df.lung<-df.lung[!is.na(df.lung$ClusteredAge),]
df.lung<-df.lung[!is.na(df.lung$group.ph.karno),]
df.lung<-df.lung[!is.na(df.lung$group.pat.karno),]

#check
sum(is.na(df.lung$ClusteredAge))+sum(is.na(df.lung$group.pat.karno))+sum(is.na(df.lung$group.ph.karno))
```

```
## [1] 0
```

Of course, nor can there  be any missing value in status or time.


```r
df.lung<-df.lung[!is.na(df.lung$status),]
df.lung<-df.lung[!is.na(df.lung$time),]
```

-------

Once we have the data, we can make a **model fitting**, either with Cox Model or Neural network, and then, we can make a **prediction**. We're going to perform the two of them.


# 1. Cox Regression model

This semiparametric model of multivariate regresion is the standard tool for estimating survival. Cox model is mainly used for analyze the importance of the covariates as prognostic factors with an easy interpretation. One of the useful facts of this model is that it also work with continuos variables in addition to categoricals, providing us a wider range of analysis.

One possible disadvantage is that the assumption of the proportionality of the hazard needs to be verified (Cox model assumes lineal influence of covariates).

We can see how  each covariate influences to the **hazard rate ** or **failure rate**, and how the survival is changed regarding to another covariate.**Hazard rate** is the probability that if somethig survives to one moment of time, it will also do for the next moment. 


**Survival modeling** can be done using Survival function $S(t)$ (survival probability at time t) or Hazard function $h(t)$ (Instantaneous hazaard rate in interval $[t,t+dt]$) give equivalent analysis. 

$h(t)=\frac{f(t)}{R(t)}=\frac{dS(t)}{dt}=h_{0}(t)·e^{ß_1x_1+ß_2x_2+...+ß_nx_n}$

- $f(t)$ Probability of density function (probability of event happening in a specific interval)

- $R(t)$ Survival function given by covariates $x_{1}, x_{2},...x{p}$, multiplied by their impact $ß{1}, ß{2},....ß{n}$ (coefficients give a idea of the importance of the covariate in the patient prognosis).

- $h_{0}(t)$ is baseline hazard (no covariates, all $x{i}=0$)


**Analyze the assigned data set and compute survival curves using the Cox regression model. Divide your data in two relevant groups, i.e., according to a relevant covariate divide the data in two groups of similar size (Age > 60, Age >=60). Analyze for these two groups the differences in survival and in the influence that the model assigns to each covariate. Use your fitted Cox model to do survival predictions for both groups and for individual patients.**


### Selection of best Cox model {#linkcovariatesselection}

Distribution of covariates that will be used: 


```r
table(df.lung$sex)
```

```
## 
##   1   2 
## 133  89
```

```r
table(df.lung$ClusteredAge)
```

```
## 
## (39,62] (62,82] 
##     102     120
```

```r
table(df.lung$group.pat.karno)
```

```
## 
##   (0,40]  (40,70] (70,100] 
##        4       74      144
```

With the `coxph()` function we're going to make a Cox model. The parameters we're going to visualize are:

- `coef` estimated coeficcient of ß (already explained).

- `exp(coef)` exponetial raised to the estimated coef.

- `se(coef)` Standard deviation of the estimation.

- `z` value of the statistical for Wald test of $ß=0$.

- `p` p-value of Wald test.


Search of best combination of covariates to explain risk function:


```r
library(FSelector)

surv.obj<-Surv(df.lung$time, df.lung$status)

exhaustiveSearch <-function(lista.vbles) {
  
  pvalue.acumulado<-list()
  
  for(nivel in 1:length(lista.vbles)) { 
    best.pvalue<-1
    for(v in lista.vbles) {
      ifelse(nivel==1, att<-v, att<-c(v, lista.acumulada))
      ec <-as.simple.formula(att,"surv.obj")
      sum<-summary(coxph(ec, data=df.lung))
      new.p<-sum$logtest["pvalue"]
      if(new.p <0.05 && new.p<best.pvalue) {
        best.pvalue<-new.p
        best.vble.nivel<-v
      }
    }
    
    lista.vbles<-lista.vbles[lista.vbles != best.vble.nivel]
    
    ifelse(nivel!=1,  pvalue.acumulado<-c(pvalue.acumulado, best.pvalue), first.pvalue<-best.pvalue)
    ifelse(nivel==1, lista.acumulada<-best.vble.nivel,lista.acumulada<-c(lista.acumulada, best.vble.nivel))
    ifelse(nivel==1, p.acum.ac<-best.pvalue, p.acum.ac<-pvalue.acumulado[nivel-1])
    if(length(lista.acumulada)!=nivel || p.acum.ac>best.pvalue[1] || length(lista.acumulada)==10) break
  }
  
df<-data.frame(lista.acumulada,unlist(list(c(first.pvalue,pvalue.acumulado))))
names(df)<-c("vble acumulada", "pvalue_acumulado")
return(df)
}

es<-exhaustiveSearch(c("age", "ClusteredAge", "sex", "group.ph.karno", "group.pat.karno", "meal.cal", "wt.loss"))

df.order<-es[order(es$pvalue_acumulado),]
df.order
```

```
##    vble acumulada pvalue_acumulado
## 2 group.pat.karno     0.0001876308
## 3         wt.loss     0.0002261348
## 5    ClusteredAge     0.0003506478
## 4             age     0.0003529887
## 6  group.ph.karno     0.0004875971
## 1             sex     0.0011895577
## 7        meal.cal     0.0054691356
```


The most important covariate for explaining the risk is the one with the groups of Karnofsky scores by the patients. Also, we can see that the clustered age is better than age itself, followed by the other score of Karnofsky (the one by the physicians), weight loss is also important, but when adding it, the pvalue gets higher.

Let's make the cox model with the top-ranked covariates and see the summary table with their respective ß and confident intervals:


```r
res.cox<-coxph(surv.obj~ClusteredAge+group.pat.karno+group.ph.karno,  data=df.lung)
summ<-summary(res.cox)
summ
```

```
## Call:
## coxph(formula = surv.obj ~ ClusteredAge + group.pat.karno + group.ph.karno, 
##     data = df.lung)
## 
##   n= 222, number of events= 161 
## 
##                             coef exp(coef) se(coef)      z Pr(>|z|)
## ClusteredAge(62,82]     -0.02958   0.97085  0.16830 -0.176    0.860
## group.pat.karno(40,70]   0.50154   1.65126  0.72602  0.691    0.490
## group.pat.karno(70,100]  0.04300   1.04394  0.73547  0.058    0.953
## group.ph.karno(40,70]    0.28094   1.32437  0.19419  1.447    0.148
## group.ph.karno(70,100]        NA        NA  0.00000     NA       NA
## 
##                         exp(coef) exp(-coef) lower .95 upper .95
## ClusteredAge(62,82]        0.9709     1.0300    0.6981     1.350
## group.pat.karno(40,70]     1.6513     0.6056    0.3979     6.852
## group.pat.karno(70,100]    1.0439     0.9579    0.2470     4.413
## group.ph.karno(40,70]      1.3244     0.7551    0.9051     1.938
## group.ph.karno(70,100]         NA         NA        NA        NA
## 
## Concordance= 0.577  (se = 0.026 )
## Rsquare= 0.055   (max possible= 0.999 )
## Likelihood ratio test= 12.51  on 4 df,   p=0.01
## Wald test            = 13.3  on 4 df,   p=0.01
## Score (logrank) test = 13.62  on 4 df,   p=0.009
```

```r
#pvalue of the model
summ$logtest[3]
```

```
##     pvalue 
## 0.01393924
```

The reference level is ClusteredAge(62,82], and regarding to it, the most influence ones are the group (40,70) of Karnofsky scores (patients unables to work but take care for most personal needs), for both, patients and physicians score. 


According to `Pr(<|z|)` colum (ratio of each regression coefficient to its standard error, a Wld statistic), `group.ph.karno(40,70]` has highly statiscally significant coeffiients

The **Hazard ratio(HZ)** is `exp(coef)` and its meaning is:

- HR=1: covariate has no effect.

- HR<1: Reduction of the Hazard probability

- HR>1: Increase of Hazard probability

So, the mentioned score groups increase the hazard probability, even more than the patients unables to take care of themselves.

Now we have the model, before using it, we have to make sure it verifys the assumption of proportional risks: *“If an individual has a risk of death at some initial time point that is twice as high as that of another individual, then at all later times the risk of death remains twice as high.” * For that purpose, we're going to use the function `cox.zph()`: the model will be accepted if $p-value \geq 0.05$.



```r
cox.zph(res.cox, transform="km", global=TRUE) 
```

```
## Warning in cor(xx, r2): the standard deviation is zero
```

```
##                             rho chisq      p
## ClusteredAge(62,82]     -0.0247 0.118 0.7310
## group.pat.karno(40,70]   0.0251 0.101 0.7501
## group.pat.karno(70,100]  0.0382 0.233 0.6291
## group.ph.karno(40,70]   -0.1313 2.906 0.0883
## group.ph.karno(70,100]       NA   NaN    NaN
## GLOBAL                       NA 6.455 0.2644
```


Global p-value is 0.2664, which means that the difference between patients is conserved in time and the **model is accepted**.


```r
plot(survfit(res.cox), xlab="Days", ylab="Portion not Rearrested")
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

Plots of the scaled Schoenfeld residuals against time. Broken lines represent the standard error band around the fit.


```r
ggc<-ggcoxzph(cox.zph(res.cox))[1:4]
library(cowplot)
```

```r
plot_grid(ggc$`1`, ggc$`4`, ggc$`3`,ggc$`2`)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


Systematic departures from a horizontal line are indicative of non-proportional hazards.



Survival curves of final choosen variables:


```r
survfit(surv.obj~ClusteredAge, df.lung, conf.type = "log-log") %>% ggsurvplot(title = "Survival by age", conf.int = T, legend.title = "Age cluster", pval=TRUE,  ggtheme = theme_bw(), risk.table=T, risk.table.col="strata", risk.table.height = 0.32)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
survfit(surv.obj~group.ph.karno, df.lung, conf.type = "log-log") %>% ggsurvplot(title = "Survival by Karnofsky scores by physics", conf.int = T, legend.title = "Score cluster", pval=TRUE,  ggtheme = theme_bw(), risk.table=T, risk.table.col="strata", risk.table.height = 0.32)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

```r
survfit(surv.obj~group.pat.karno, df.lung, conf.type = "log-log") %>% ggsurvplot(title = "Survival by Karnofsky scores by patients", conf.int = T, legend.title = "Score cluster", pval=TRUE,  ggtheme = theme_bw(), risk.table=T, risk.table.col="strata", risk.table.height = 0.32)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-18-3.png)<!-- -->
 
 The biggest difference between levels is found in the covariate *Karfnosky scores by physicians* ($p-value=0.0073$), followed by *Karfnosky scores by physicians* (the plot is distort by the level (0,40] which has a wider confidence interval due to the lack of samples in this group, only 4). Clustered age is the worst one, let's pick another level as the reference one for the cox model to see the ß value:
 

```r
res.cox2<-coxph(surv.obj~group.pat.karno+group.ph.karno+ClusteredAge,  data=df.lung)
summ2<-summary(res.cox2)
summ2
```

```
## Call:
## coxph(formula = surv.obj ~ group.pat.karno + group.ph.karno + 
##     ClusteredAge, data = df.lung)
## 
##   n= 222, number of events= 161 
## 
##                             coef exp(coef) se(coef)      z Pr(>|z|)
## group.pat.karno(40,70]   0.50154   1.65126  0.72602  0.691    0.490
## group.pat.karno(70,100]  0.04300   1.04394  0.73547  0.058    0.953
## group.ph.karno(40,70]    0.28094   1.32437  0.19419  1.447    0.148
## group.ph.karno(70,100]        NA        NA  0.00000     NA       NA
## ClusteredAge(62,82]     -0.02958   0.97085  0.16830 -0.176    0.860
## 
##                         exp(coef) exp(-coef) lower .95 upper .95
## group.pat.karno(40,70]     1.6513     0.6056    0.3979     6.852
## group.pat.karno(70,100]    1.0439     0.9579    0.2470     4.413
## group.ph.karno(40,70]      1.3244     0.7551    0.9051     1.938
## group.ph.karno(70,100]         NA         NA        NA        NA
## ClusteredAge(62,82]        0.9709     1.0300    0.6981     1.350
## 
## Concordance= 0.577  (se = 0.026 )
## Rsquare= 0.055   (max possible= 0.999 )
## Likelihood ratio test= 12.51  on 4 df,   p=0.01
## Wald test            = 13.3  on 4 df,   p=0.01
## Score (logrank) test = 13.62  on 4 df,   p=0.009
```

In effect, it's the poorest of the influences.

We can also see the percentage of events occured in each selected group:



```r
with(df.lung, prop.table(table(ClusteredAge, status)))
```

```
##             status
## ClusteredAge         1         2
##      (39,62] 0.1396396 0.3198198
##      (62,82] 0.1351351 0.4054054
```

```r
with(df.lung, prop.table(table(group.pat.karno, status)))
```

```
##                status
## group.pat.karno           1           2
##        (0,40]   0.009009009 0.009009009
##        (40,70]  0.058558559 0.274774775
##        (70,100] 0.207207207 0.441441441
```

```r
with(df.lung, prop.table(table(group.ph.karno, status)))
```

```
##               status
## group.ph.karno          1          2
##       (0,40]   0.00000000 0.00000000
##       (40,70]  0.03153153 0.21621622
##       (70,100] 0.24324324 0.50900901
```

The groups with more deaths are (70,100] in physicians and patients score, followed by (62,82] age group.

----

### Cox model as a predictor 



Type of predictions values:

- `"risk"`: Predict the risk score ($exp(lp)$). 

- `"expected"`: Predict the expected number of events given the covariates and follow-up time.

- Survival probability = $exp(-expected)$

- `"terms"`: Predict the terms of the linear predictor


We're going to use and compare just the survival probability.

#### With 3 covariates

First we're going to train and predict the survival probability in time using the **same data** (all the dataframe). This is not really useful at all, because if the model is trained with the same data that is going to predict will learn it, and besides the result is not "real", we won't see the ability to generalize (what really matters in prediction models) but we're going to see how it works in this trial. `coxph()` will train the model with the selected variables, and `predict(cox.model, type=expected)` will do the predictions from the trained Cox model.


```r
#Using all the data for training the cox model 
cox.model<-coxph(Surv(time,status)~ClusteredAge+group.pat.karno+group.ph.karno,  data=df.lung)
predictions<-as.vector(sort(exp(-predict(cox.model, type=c("expected"))),decreasing=T))
survfit.values<-sort(survfit(formula=surv.obj~ClusteredAge+group.pat.karno+group.ph.karno, data=df.lung)$surv, decreasing=T)
all.times<-sort(df.lung$time)
```






```r
plot(all.times, predictions, ylim=c(0,1),col="red",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(all.times,survfit.values, ylim=c(0,1),col="blue",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("real","predicted"),lty=c(1,1),col=c(2,4))
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

As could not be otherwise (for the reasons explained above) the predicted and the `survfut`data are very similar, but couldn't be exactly the same due to the low amount of data the model has been trained with.

The correct way to do train models and do predictions is making a split in the data, and verify the quality of predictions with new data the model has never seen before. This isn't correct at all (like a Cross validation or other similar method would be) cause with only one split, the results can be due to this particular separation, but just for a trial.

As before, `coxph` will train the model just with the train data, but this time, `survfit` with the trained Cox model has the formula parameter will make the predictions over the test data.


```r
train<-df.lung[1:180,]
test<-df.lung[181:222,]

cox.model<-coxph(Surv(time,status)~ClusteredAge+group.pat.karno+group.ph.karno,  data=train)
#Computes the predicted survivor function for a Cox proportional hazards model
predictions<-survfit(formula=cox.model, newdata=test) 
```

This way neither reach our purpose at all because it gives the predicted probability for each patient (row of test data) for each time:


```r
head(summary(predictions)$surv) 
```

```
##            185       186       187       188       189       190       191
## [1,] 0.9956136 0.9956136 0.9946643 0.9956136 0.9934065 0.9956136 0.9957171
## [2,] 0.9868397 0.9868397 0.9840070 0.9868397 0.9802620 0.9868397 0.9871491
## [3,] 0.9824499 0.9824499 0.9786826 0.9824499 0.9737076 0.9824499 0.9828616
## [4,] 0.9780430 0.9780430 0.9733425 0.9780430 0.9671423 0.9780430 0.9785569
## [5,] 0.9736185 0.9736185 0.9679865 0.9736185 0.9605660 0.9736185 0.9742346
## [6,] 0.9691839 0.9691839 0.9626235 0.9691839 0.9539897 0.9691839 0.9699020
##            192       193       194       195       196       197       198
## [1,] 0.9915898 0.9956136 0.9957171 0.9957171 0.9957171 0.9957171 0.9957171
## [2,] 0.9748697 0.9868397 0.9871491 0.9871491 0.9871491 0.9871491 0.9871491
## [3,] 0.9665555 0.9824499 0.9828616 0.9828616 0.9828616 0.9828616 0.9828616
## [4,] 0.9582430 0.9780430 0.9785569 0.9785569 0.9785569 0.9785569 0.9785569
## [5,] 0.9499321 0.9736185 0.9742346 0.9742346 0.9742346 0.9742346 0.9742346
## [6,] 0.9416370 0.9691839 0.9699020 0.9699020 0.9699020 0.9699020 0.9699020
##            199       200       201       202       203       204       205
## [1,] 0.9956136 0.9932472 0.9957171 0.9917880 0.9957171 0.9956136 0.9957171
## [2,] 0.9868397 0.9797885 0.9871491 0.9754570 0.9871491 0.9868397 0.9871491
## [3,] 0.9824499 0.9730790 0.9828616 0.9673338 0.9828616 0.9824499 0.9828616
## [4,] 0.9780430 0.9663595 0.9785569 0.9592106 0.9785569 0.9780430 0.9785569
## [5,] 0.9736185 0.9596298 0.9742346 0.9510874 0.9742346 0.9736185 0.9742346
## [6,] 0.9691839 0.9529013 0.9699020 0.9429779 0.9699020 0.9691839 0.9699020
##            207       208       209       210       211       212       213
## [1,] 0.9956136 0.9915898 0.9934065 0.9956136 0.9956136 0.9934065 0.9934065
## [2,] 0.9868397 0.9748697 0.9802620 0.9868397 0.9868397 0.9802620 0.9802620
## [3,] 0.9824499 0.9665555 0.9737076 0.9824499 0.9824499 0.9737076 0.9737076
## [4,] 0.9780430 0.9582430 0.9671423 0.9780430 0.9780430 0.9671423 0.9671423
## [5,] 0.9736185 0.9499321 0.9605660 0.9736185 0.9736185 0.9605660 0.9605660
## [6,] 0.9691839 0.9416370 0.9539897 0.9691839 0.9691839 0.9539897 0.9539897
##            214       215       216       217       218       219       220
## [1,] 0.9946643 0.9957171 0.9934065 0.9956136 0.9934065 0.9933603 0.9956136
## [2,] 0.9840070 0.9871491 0.9802620 0.9868397 0.9802620 0.9801248 0.9868397
## [3,] 0.9786826 0.9828616 0.9737076 0.9824499 0.9737076 0.9735255 0.9824499
## [4,] 0.9733425 0.9785569 0.9671423 0.9780430 0.9671423 0.9669155 0.9780430
## [5,] 0.9679865 0.9742346 0.9605660 0.9736185 0.9605660 0.9602947 0.9736185
## [6,] 0.9626235 0.9699020 0.9539897 0.9691839 0.9539897 0.9536743 0.9691839
##            221       222       223       224       226       227       228
## [1,] 0.9957171 0.9957171 0.9957171 0.9934065 0.9917880 0.9957171 0.9956136
## [2,] 0.9871491 0.9871491 0.9871491 0.9802620 0.9754570 0.9871491 0.9868397
## [3,] 0.9828616 0.9828616 0.9828616 0.9737076 0.9673338 0.9828616 0.9824499
## [4,] 0.9785569 0.9785569 0.9785569 0.9671423 0.9592106 0.9785569 0.9780430
## [5,] 0.9742346 0.9742346 0.9742346 0.9605660 0.9510874 0.9742346 0.9736185
## [6,] 0.9699020 0.9699020 0.9699020 0.9539897 0.9429779 0.9699020 0.9691839
```

We can also see the survival probability value for all patients for a certain time, for example t=50 and t=1000.

```r
summary(predictions, times=50)
```

```
## Call: survfit(formula = cox.model, newdata = test)
## 
##  time n.risk n.event survival1 survival2 survival3 survival4 survival5
##    50    171       9      0.96      0.96     0.952      0.96     0.941
##  survival6 survival7 survival8 survival9 survival10 survival11 survival12
##       0.96     0.961     0.925      0.96      0.961      0.961      0.961
##  survival13 survival14 survival15 survival16 survival17 survival18
##       0.961      0.961       0.96      0.939      0.961      0.927
##  survival19 survival20 survival21 survival22 survival23 survival24
##       0.961       0.96      0.961       0.96      0.925      0.941
##  survival25 survival26 survival27 survival28 survival29 survival30
##        0.96       0.96      0.941      0.941      0.952      0.961
##  survival31 survival32 survival33 survival34 survival35 survival36
##       0.941       0.96      0.941       0.94       0.96      0.961
##  survival37 survival38 survival39 survival40 survival41 survival42
##       0.961      0.961      0.941      0.927      0.961       0.96
```

```r
summary(predictions, times=1000)
```

```
## Call: survfit(formula = cox.model, newdata = test)
## 
##  time n.risk n.event survival1 survival2 survival3 survival4 survival5
##  1000      2     147    0.0871    0.0871    0.0513    0.0871    0.0254
##  survival6 survival7 survival8 survival9 survival10 survival11 survival12
##     0.0871    0.0923    0.0092    0.0871     0.0923     0.0923     0.0923
##  survival13 survival14 survival15 survival16 survival17 survival18
##      0.0923     0.0923     0.0871     0.0232     0.0923     0.0103
##  survival19 survival20 survival21 survival22 survival23 survival24
##      0.0923     0.0871     0.0923     0.0871     0.0092     0.0254
##  survival25 survival26 survival27 survival28 survival29 survival30
##      0.0871     0.0871     0.0254     0.0254     0.0513     0.0923
##  survival31 survival32 survival33 survival34 survival35 survival36
##      0.0254     0.0871     0.0254     0.0248     0.0871     0.0923
##  survival37 survival38 survival39 survival40 survival41 survival42
##      0.0923     0.0923     0.0254     0.0103     0.0923     0.0871
```



#### Just with covariate of age division

These previous results are a lot of data to visualize or clearlu understand anything. For simplicity, we're going to keep only one of the variables for the Cox model, and make a different dataframe for each level of it.


```r
cox.age<-coxph(Surv(time,status)~ClusteredAge,  data=df.lung)
summ<-summary(cox.age)
summ$logtest[3]
```

```
##    pvalue 
## 0.4211291
```

```r
cox.zph(cox.age, transform="km", global=TRUE)
```

```
##                         rho chisq     p
## ClusteredAge(62,82] -0.0591 0.554 0.457
```

The pvalue of the model using just this covarite is a little bit worse comparing with the previous set of variables (0.42 compared to 0.01). The validity test gives the same result.

We're going to make a Cox model and a prediction for each level for all the patients for the selected time. Also, we will obtain the `survfit` (real) value (the mean of all the patients in the same time). 

### Predicted curve vs Real curve bu clusters


```r
library(pec)
library(riskRegression)
```


```r
set.seed(123)
prediction.and.reference<-function(all.data,time, level) {
  
  #Make a dataframe for each level 
  df<-dplyr::filter(all.data, all.data$ClusteredAge==level)
  #Split the data in train and test
  train_ind <- sample(seq_len(nrow(df)), size = floor(0.8*nrow(df)))
  train<-df[train_ind, ]
  test<-df[-train_ind, ]
  #Create cox model
  cox.model<-coxph(Surv(time,status)~1,  data=train,x=T)
  #Predict survival data for time using cox model
  predictions<-predictSurvProb(cox.model,newdata=test, times=time)[,1]
  #Obtain "real" values from survfit function for time
  reference.values<-survfit(Surv(time,status)~1, data=test)
  summ.time<-summary(reference.values, times=time)$surv
  return(list(predictions,summ.time, mean(predictions)))
}
```


```r
levels(df.lung$ClusteredAge)
```

```
## [1] "(39,62]" "(62,82]"
```

```r
#CLUSTER1 
c1<-prediction.and.reference(df.lung,95,"(39,62]")
#Predictions for all patients in time 95
c1[[1]]
```

```
##  [1] 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813
##  [8] 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813
## [15] 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813 0.9263813
```

```r
#Mean of predictions for time 95
c1[[3]]
```

```
## [1] 0.9263813
```

```r
#Reference value for time 95 (Survfit value)
c1.reference<-c1[[2]]
#CLUSTER2
c2<-prediction.and.reference(df.lung,95,"(62,82]")
c2[[1]]
```

```
##  [1] 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359
##  [8] 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359
## [15] 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359 0.8238359
## [22] 0.8238359 0.8238359 0.8238359
```

```r
c2[[3]]
```

```
## [1] 0.8238359
```

```r
c2[[2]]
```

```
## [1] 0.8333333
```




```r
df.lung.cluster1<-dplyr::filter(df.lung, df.lung$ClusteredAge=="(39,62]")
df.lung.cluster2<-dplyr::filter(df.lung, df.lung$ClusteredAge=="(62,82]")
tiempos1<-sort(df.lung.cluster1$time)
tiempos2<-sort(df.lung.cluster2$time)
c1.predicted.mean<-sapply(tiempos1, function(t) prediction.and.reference(df.lung.cluster1,t,"(39,62]")[[3]])
c2.predicted.mean<-sapply(tiempos2, function(t) prediction.and.reference(df.lung.cluster2,t,"(62,82]")[[3]])
c1.reference.mean<-sapply(tiempos1, function(t) prediction.and.reference(df.lung.cluster1,t,"(39,62]")[[2]])
c2.reference.mean<-sapply(tiempos2, function(t) prediction.and.reference(df.lung.cluster2,t,"(62,82]")[[2]])
```










```r
c1.predicted.mean
```

```
##   [1] 0.98773022 0.97546043 0.98773022 0.96319065 0.95092088 0.93865110
##   [7] 0.91411156 0.91411156 0.88957202 0.90184179 0.88957202 0.87730226
##  [13] 0.82822324 0.85276274 0.86503250 0.85276274 0.84049299 0.84049299
##  [19] 0.80368375 0.80368375 0.81595349 0.77914427 0.76687453 0.77914427
##  [25] 0.80368375 0.75460481 0.72942523 0.76567754 0.75342696 0.72884859
##  [31] 0.66606408 0.70387631 0.70277254 0.67754002 0.68827518 0.66760718
##  [37] 0.63675132 0.70116653 0.66532360 0.71299057 0.62473648 0.67250543
##  [43] 0.64960167 0.63410975 0.65756650 0.59715754 0.59596030 0.59537478
##  [49] 0.64618761 0.57253677 0.55431223 0.61661907 0.53132417 0.54162049
##  [55] 0.56792340 0.58464666 0.58712114 0.55436679 0.50562631 0.53784746
##  [61] 0.53729330 0.52387303 0.53172758 0.52576031 0.49902924 0.52279750
##  [67] 0.49399263 0.44507903 0.40731955 0.47729161 0.46500093 0.42919134
##  [73] 0.39219950 0.38579254 0.35720927 0.43520078 0.35834363 0.35394824
##  [79] 0.38851794 0.36750004 0.34430487 0.35451492 0.38271247 0.33890140
##  [85] 0.32933492 0.35269695 0.28505103 0.29905785 0.27757173 0.26691104
##  [91] 0.27390500 0.23724728 0.21419344 0.17633869 0.19721446 0.11745073
##  [97] 0.11868718 0.12067091 0.06412818 0.07779680 0.02414477 0.05177099
```

```r
c1.reference.mean
```

```
##   [1] 1.00000000 0.95238095 1.00000000 0.90476190 0.95238095 1.00000000
##   [7] 0.95238095 0.90476190 0.95238095 0.85714286 0.80952381 0.80952381
##  [13] 0.90476190 0.85714286 0.80952381 0.85714286 0.76190476 0.71428571
##  [19] 0.85714286 0.76190476 0.71428571 0.80952381 0.80952381 0.80952381
##  [25] 0.90476190 0.76190476 0.80952381 0.76190476 0.71428571 0.57142857
##  [31] 0.66666667 0.61904762 0.76190476 0.65476190 0.90476190 0.61904762
##  [37] 0.60952381 0.61538462 0.76190476 0.56122449 0.75294118 0.70158730
##  [43] 0.70833333 0.80423280 0.65641026 0.66666667 0.65476190 0.75000000
##  [49] 0.58543417 0.57142857 0.59940060 0.50595238 0.70370370 0.57142857
##  [55] 0.53968254 0.47165533 0.60084034 0.65934066 0.44698413 0.54421769
##  [61] 0.68399235 0.47142857 0.59369202 0.65384615 0.61904762 0.42631173
##  [67] 0.54545455 0.59728507 0.24928775 0.41904762 0.37037037 0.35374150
##  [73] 0.36940837 0.31730769 0.27428571 0.57880952 0.30919312 0.38480038
##  [79] 0.34632035 0.39899554 0.27500000 0.25972222 0.31250000 0.27777778
##  [85] 0.31250000 0.20634921 0.07857143 0.12155745 0.32980600 0.26069519
##  [91] 0.17000000 0.41958042 0.23376623 0.11000000 0.34151786 0.15462896
##  [97] 0.11446886 0.09308339 0.10000000 0.08000000 0.08000000 0.07000000
```

```r
plot(tiempos1, c1.predicted.mean, ylim=c(0,1),col="red",main="Survival curves for age cluster (39,62]",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(tiempos1,c1.reference.mean,ylim=c(0,1),col="blue",main="Survival curves for age cluster (39,62]",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("real","predicted"),lty=c(1,1),col=c(2,4))
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

```r
c2.predicted.mean
```

```
##   [1] 1.00000000 0.96891220 0.96891220 0.96891220 0.95854960 0.93782441
##   [7] 0.92746181 0.92746181 0.90673663 0.92746181 0.91709922 0.91709922
##  [13] 0.88601145 0.89637404 0.86528627 0.87564886 0.85492369 0.85492369
##  [19] 0.84456110 0.84456110 0.85492369 0.81334136 0.82383594 0.80270960
##  [25] 0.81295182 0.82294267 0.78198452 0.76018258 0.79157002 0.75049569
##  [31] 0.79120937 0.77089263 0.73976228 0.74093537 0.74883900 0.71742259
##  [37] 0.75115500 0.73740917 0.69536874 0.71844365 0.69507663 0.66396114
##  [43] 0.67436289 0.68620770 0.65379140 0.65173739 0.68317017 0.62303758
##  [49] 0.66265840 0.66138930 0.61827390 0.61869400 0.64287814 0.64042115
##  [55] 0.62034208 0.59603939 0.64193458 0.58195847 0.58634282 0.61629746
##  [61] 0.62649853 0.56043941 0.60963438 0.59277605 0.60240096 0.56604249
##  [67] 0.58213814 0.56507967 0.57626040 0.51122568 0.53806528 0.49144009
##  [73] 0.52514252 0.51808311 0.48995812 0.45968221 0.48853985 0.47021945
##  [79] 0.40465965 0.44266363 0.42412297 0.43164575 0.38988655 0.37599069
##  [85] 0.41404212 0.39796605 0.35165728 0.34896900 0.39388735 0.35600616
##  [91] 0.34926284 0.33396679 0.30555100 0.28276191 0.32541304 0.27923935
##  [97] 0.32554395 0.31973469 0.29661980 0.30209431 0.26354476 0.23066211
## [103] 0.21694565 0.25025164 0.24042212 0.19468438 0.20072234 0.17933668
## [109] 0.16642439 0.11969943 0.12198494 0.12399950 0.12004201 0.09826411
## [115] 0.07607305 0.04194973 0.08265342 0.06943613 0.06730102 0.08166130
```

```r
c2.reference.mean
```

```
##   [1] 1.00000000 0.95833333 0.95833333 1.00000000 0.91666667 0.95833333
##   [7] 1.00000000 0.91666667 1.00000000 0.87500000 0.95833333 0.87500000
##  [13] 0.75000000 0.87500000 0.83333333 0.91666667 0.91666667 0.87500000
##  [19] 0.70833333 0.75000000 0.75000000 0.83333333 0.79166667 0.91666667
##  [25] 0.91666667 0.91666667 0.78947368 0.66666667 0.70833333 0.70833333
##  [31] 0.74561404 0.70833333 0.70833333 0.70833333 0.91477273 0.83333333
##  [37] 0.75000000 0.58333333 0.70000000 0.70175439 0.75000000 0.66176471
##  [43] 0.74375000 0.74375000 0.65789474 0.54166667 0.66666667 0.70588235
##  [49] 0.74375000 0.70175439 0.58333333 0.69568452 0.73569024 0.65625000
##  [55] 0.54700855 0.48379630 0.56919643 0.66176471 0.48951049 0.62222222
##  [61] 0.56538462 0.52083333 0.64615385 0.65333333 0.56286550 0.55654762
##  [67] 0.55151515 0.47088068 0.45454545 0.35108025 0.46897547 0.42849794
##  [73] 0.53199405 0.45000000 0.45117187 0.48375154 0.42797203 0.49845679
##  [79] 0.49242424 0.55416667 0.59375000 0.48897059 0.58333333 0.39111111
##  [85] 0.53472222 0.49218750 0.41894531 0.34090909 0.36060606 0.19308894
##  [91] 0.24305556 0.19191919 0.30864198 0.32096531 0.23034188 0.50854701
##  [97] 0.37799640 0.52500000 0.23123605 0.33124770 0.29501748 0.28125000
## [103] 0.20409689 0.43928683 0.13773148 0.31100478 0.27562500 0.08000000
## [109] 0.05989859 0.38287450 0.06324405 0.15028831 0.05252101 0.06111111
## [115] 0.14000000 0.05902778 0.20000000 0.15000000 0.12000000 0.08000000
```

```r
plot(tiempos2, c2.predicted.mean, ylim=c(0,1),col="red",main="Survival curves for age cluster (62,82]",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(tiempos2,c2.reference.mean,ylim=c(0,1),col="blue",main="Survival curves for age cluster (62,82]",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("real","predicted"),lty=c(1,1),col=c(2,4))
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-36-2.png)<!-- -->

The curves aren't too smooth due to the insuficient amount of data we have, we only used one level in each dataframe (it's like more or less the half of the original data, which was already few). As expected, the predicted one has even more peaks. We have to consider that in these conditions, a difference of just a few patients is very notable in the curve. Regardless of this, I would see the prediction is good.




-----


#2. NNET


Neural networls are computer models that matches the functionality of the brain in a very simple manner. For that, they use neurons (elementary processing unitof the nervous system). 

They way the work very briefly:

- **Biological neural network**: Signals are received at the dendrites from other neurons through the synapses, and they are combined in the soma to yield an electric potential. Signals are transmitted along the axon, and they depend on the potential.


- **Artificial neural network**: Used to solve artificial intelligence problems. The connections of networks are modelled as weigths. The inputs are modified by these weights and summed. In this process the neurons "learn"" (update or modify the connections weights in response to stimuli being presented at the input layer and optionally the output layer).


### Data replication 

Algorithms of statistical learning can’t learn well from a few examples. Neural networks model the underlying distribution of the dataset that they are being trained on. The dataset itself, however big it is, is only a tiny fraction of the general population, and in nature there's always fluctuations, not only a strict pattern. Our goal is to obtain the pattern of all population with a sample N. The bigger N, the "less the differences", the nearest to find  a pattern.

Knowing this all, our previous dataset was too small and insufficient, for having more samples using the same data and values, for each patient we're going to replicate the colum values for a wide range of times, mantaining the state in each time.


```r
last.time<-max(df.lung$time)
replicate.df<-data.frame()
 
for (patient in 1:nrow(df.lung)) {
  cols.to.save<-c("inst", "age", "sex", "ph.ecog", "ph.karno", "pat.karno", "meal.cal", "wt.loss", "ClusteredAge", "group.ph.karno", "group.pat.karno")
  cols<-sapply(cols.to.save, function(c) df.lung[,c][patient])
  status<-df.lung$status[patient]
  time<-df.lung$time[patient]
  #status alive
 for(row in seq(1,time,by=5)) {
    new.row<-c(cols[1], row, 1, cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cols[11])
    replicate.df<-rbind(replicate.df, new.row) 
  } 
for(row in seq(time,last.time,by=5)) {
  if(status==1) {
    new.row<-c(cols[1], row, 1, cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cols[11])
  }
  else if(status==2) {
    new.row<-c(cols[1], row, 2, cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols[8], cols[9], cols[10], cols[11])
  }
  replicate.df<-rbind(replicate.df, new.row) 
}
}
 names(replicate.df)<-colnames(df.lung)
```







```r
dim(df.lung)
```

```
## [1] 222  13
```

```r
dim(replicate.df)
```

```
## [1] 45598    13
```

Our previous dataset had 222 samples. This new one has 14033, not too large yet but enough for our purpose without having too slow execution.



### Fit the NNET model 


Before using the model, we have to find and define the **most suitable parameters and architecture** according to our data.

The complexity of our **architecture** should be choosen according to the complexity of our data. Since has we have just seen it is pretty simple, we're going to use a Multilayer Perceptron with just **one middle hidden layer**: Feedforward network. The network will work with Forward activity, receiving inputs (train data) and giving outputs (predictions). We will later see in the network plot.


Knowing the number of layers, we still have to choose the best **parameters** of **size** (number of networks) and **decay** yo prevent overfitting effects (if there are so many neurons the model learn the train data so good and lost the ability to generalize).

#### Train and test data

For validating the results of predictions and see the ability to generalize of our model (our ultimate goal). We need to save apart some data that the model has never seen before during the training (the set test). 

First we split the data, using 80% of randomly choosen samples for train and the remaining 20% will be test. Then, we're going to make the Class (the one to predict) with the status values. The NA values will be imputed with the mean.



```r
library(nnet)


train_ind <- sample(seq_len(nrow(replicate.df)), size = floor(0.8*nrow(replicate.df)))
nnet.train<-replicate.df[train_ind, ]
nnet.test<-replicate.df[-train_ind, ]

nnet.train$ClusteredAge<-factor(nnet.train$ClusteredAge, labels=c("(39,62]","(62,82]"))
nnet.test$ClusteredAge<-factor(nnet.test$ClusteredAge, labels=c("(39,62]","(62,82]"))
  
#Create class colum and remove status
nnet.train$Class <- class.ind(nnet.train$status)
nnet.train<-nnet.train[ , !(names(nnet.train) %in% c("status"))]
nnet.test$Class <- class.ind(nnet.test$status)
nnet.test<-nnet.test[ , !(names(nnet.test) %in% c("status"))]

#Impute missing values with mean
for(col in 1:ncol(nnet.train)) nnet.train[is.na(nnet.train[,col]), col] <- mean(nnet.train[,col], na.rm = TRUE)
for(col in 1:ncol(nnet.test)) nnet.test[is.na(nnet.test[,col]), col] <- mean(nnet.test[,col], na.rm = TRUE)

replicate.df<-rbind(nnet.train, nnet.test)
```

#### Best parameters

Before fitting the model to predict the survival data we need to explore and discover which parameters are the most suitable ones for NNET algorithm using our data. For this purpose, we will check the AUC performance value of different combinations of:

- `decay:` weight decay: Uses as penalty the sum of squares of the weights $w_{ij}$. The use of weight decay seems both to help the optimization process and to avoid over-ﬁtting. 

- `size:` Number of units (neurons) in the hidden layer.

To rapidly find the best combination of size and decay in our data, we're going to see the AUC result of all the posible combinations of two vectors of values, having a wide range of posibilities.

In this first step, we are using `type="class"`. With this adjustement, the outputs will be $0$ and $1$ (deceased or living for that time), that we're going to compare directly with the Class of test (status) to see how good it works before predicting the survival probabilities.






```r
library(pROC)

auc.for.parameters<-function(train, test, sizes, decays) {
  all.aucs<-vector()
  l.sizes<-list()
  l.decays<-list()
  for(d in decays) {
    for(s in sizes) {
      model<-nnet(Class~. ,size=s, data=train, softmax = TRUE, entropy = TRUE, decay = d)
      predictions<-predict(model, test, type = "class")
      auc<-roc(test$Class[,"1"], as.numeric(predictions), smooth=FALSE, auc=TRUE)$auc
      all.aucs<-c(all.aucs, auc)
      l.sizes<-c(l.sizes, s)
      l.decays<-c(l.decays, d)
    }
  }
  df.combinations<-data.frame(unlist(l.decays), unlist(l.sizes), all.aucs)
  names(df.combinations)<-c("decay", "size", "auc")
 return(df.combinations)
}

sizes=c(1,10,15,20,25,30,50)
decays=c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
combinations<-auc.for.parameters(nnet.train, nnet.test, sizes, decays)
```




 
Visualizing only the first 10 top-ranked values we have:


```r
combinations.sortedby.auc<-combinations[order(combinations$auc, decreasing = TRUE),][1:10,]
combinations.sortedby.auc
```

```
##    decay size       auc
## 12  0.05   25 0.7745288
## 2   0.01   10 0.7677460
## 24  0.20   15 0.7658349
## 5   0.01   25 0.7566555
## 35  0.30   50 0.7526153
## 17  0.10   15 0.7495223
## 27  0.20   30 0.7467362
## 30  0.30   10 0.7461851
## 3   0.01   15 0.7443197
## 18  0.10   20 0.7441337
```


```r
library("scatterplot3d")
s3d<-scatterplot3d(x=combinations.sortedby.auc,  xlab = "#neurons", ylab = "decay", color="steelblue", pch = 16, box=FALSE, grid =TRUE,type="h")
s3d.coords <-s3d$xyz.convert(combinations.sortedby.auc)
text(s3d.coords$x, s3d.coords$y, labels=combinations.sortedby.auc[,3], cex=0.7, pos=3)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-48-1.png)<!-- -->


### Survival probabilities prediction

Best parameters: **decay=0.05**, **size=25**. The other good values are proportional, smaller number and neurons and decay. If there are less neurons learning, ther's less overfitting control needed.

Now we now that the model works well with that parameters, we're going to use `type="raw"` to predict the survival probabilities. 

As formula, we will use the same as in the Cox model for a more accurate comparision (`Class~ClusteredAge`). 


```r
model<-nnet(Class~.,size=25,data=nnet.train, softmax = TRUE, entropy = TRUE, decay = 0.05)
predictions<-predict(model, nnet.test, type = "raw")
```





The result of predictons will be a matrix of the same size of Class, having for each row (sample) the probability of be one class or another (1 or 2):



```r
predictions[1:20,]
```

```
##            1          2
## 19 0.9244296 0.07557040
## 23 0.9241637 0.07583626
## 25 0.9240133 0.07598671
## 26 0.9239332 0.07606676
## 28 0.9237628 0.07623716
## 33 0.9232696 0.07673039
## 36 0.9229208 0.07707916
## 45 0.7194001 0.28059991
## 50 0.7166358 0.28336423
## 51 0.7160142 0.28398580
## 53 0.7132553 0.28674469
## 57 0.5344698 0.46553023
## 59 0.5324537 0.46754630
## 63 0.5291381 0.47086186
## 64 0.5279504 0.47204962
## 78 0.5061936 0.49380644
## 81 0.5000975 0.49990252
## 82 0.4979396 0.50206044
## 87 0.4861686 0.51383143
## 92 0.4727148 0.52728516
```

### Predicted curve vs Real curve


For making the curve with survfit in order to compare it with te "real" one, we're going to obtain the "predicted class" choosing the one with the bigger probability.


```r
predicted.status<-vector()
for(row in 1:nrow(predictions)) {
  censored.prob<-predictions[row,1]
  dead.prob<-predictions[row,2]
  predicted.status[row]<-ifelse(censored.prob>dead.prob, 0,1)
}

nnet.test$predicted.status<-predicted.status
```



Using the same formula of Cox Model:


```r
surv.predicted<-survfit(Surv(time, predicted.status)~ClusteredAge, data=nnet.test)
surv.real<-survfit(Surv(time, Class[,1])~ClusteredAge, data=nnet.test)

splots<-list()
splots[[1]]<-ggsurvplot(surv.real, conf.int = TRUE, censor= TRUE, cex.axis=3, cex.lab=3.0, main="Real Survival curve Age grouping", pval = TRUE, palette = "Dark2")
splots[[2]]<-ggsurvplot(surv.predicted, conf.int = TRUE, censor= TRUE, cex.axis=3, cex.lab=3.0, main="Predicted Survival curve Age grouping", pval = TRUE, palette = "Dark2")
arrange_ggsurvplots(splots, print = TRUE, ncol =2, nrow = 1, risk.table.height = 1)
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

The p-value of the predicted curve is as goo as the real one. The shape is normal, but in order to see the real difference between predicted and real we will see them separated by the clusters:


### Predicted curve vs Real curve by clusters



```r
replicate.df.cluster1<-dplyr::filter(replicate.df, replicate.df$ClusteredAge=="(39,62]")
replicate.df.cluster2<-dplyr::filter(replicate.df, replicate.df$ClusteredAge=="(62,82]")
tiempos1<-sort(replicate.df.cluster1$time)
tiempos2<-sort(replicate.df.cluster2$time)

cluster.predictions<-function(all.data, level) {
  
  #Make a dataframe for each level 
  df<-dplyr::filter(all.data, all.data$ClusteredAge==level)
  #Split the data in train and test
  train_ind <- sample(seq_len(nrow(df)), size = floor(0.8*nrow(df)))
  train<-df[train_ind, ]
  test<-df[-train_ind, ]
  #Create model
  model<-nnet(Class~.,size=25,data=train, softmax = TRUE, entropy = TRUE, decay = 0.05)
  #Predict survival probabilities
  predictions<-predict(model, test, type = "raw")
  predicted.status<-vector()
  
  for(row in 1:nrow(predictions)) {
    censored.prob<-predictions[row,1]
    dead.prob<-predictions[row,2]
    predicted.status[row]<-ifelse(censored.prob>dead.prob, 0,1)
  }
  test$predicted.status<-predicted.status
  surv.predicted<-survfit(Surv(time, predicted.status)~ClusteredAge, data=test)
  surv.real<-survfit(Surv(time, Class[,1])~ClusteredAge, data=test)
  t<-test$time
  return(list(surv.predicted, surv.real, t))
}
debug(cluster.predictions)
c1<-cluster.predictions(replicate.df, "(39,62]")
c2<-cluster.predictions(replicate.df,  "(62,82]")
```








```r
par(mfrow=c(1,2))

plot(c1[[2]], col=2)
par(new=TRUE)
plot(c1[[1]], col=4)
legend("topright",c("(39,62] real","(39,62] predicted"),lty=c(1,1),col=c(2,4))

plot(c2[[2]], col=2)
par(new=TRUE)
plot(c2[[1]], col=4)
legend("topright",c("(62,82] real","(62,82] predicted"),lty=c(1,1),col=c(2,4))
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-56-1.png)<!-- -->

The predicted curve show more probabolity of survival at begining than real. For the second level of age , the predicted curve shows less survival at the end, probably because of the split of the data. 


#3. Comparasions


We have already seen all the comparations by clusters for each model and dataset.

Using both levels:

**Real survival curve with lung dataset**


```r
#Original lung dataframe
surv.real.lung<-survfit(Surv(time, status)~ClusteredAge, data=df.lung)
ggsurvplot(surv.real.lung, conf.int = TRUE, censor= TRUE, cex.axis=3, cex.lab=3.0, pval = TRUE, main="Real curve lung dataset")
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-57-1.png)<!-- -->


**Real survival curve with replicated dataset**


```r
#Replicated dataframe
surv.real.replicate<-survfit(Surv(time, Class[,1])~ClusteredAge, data=replicate.df)
ggsurvplot(surv.real.replicate, conf.int = TRUE, censor= TRUE, cex.axis=3, cex.lab=3.0, pval = TRUE, main="Real curve replicated dataset")
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-58-1.png)<!-- -->


```r
dim(df.lung)
```

```
## [1] 222  13
```

```r
dim(replicate.df)
```

```
## [1] 45598    13
```

Regarding to the real curves: 

The one using the replicated dataset has a much better $p-value$ and the confident intervals are too small that are almost negligible in the curve. That is because using much more data give more reliable results. Also, the end of the curves are sighlty different, being the survival considerably  higher in the curve that uses the replicated dataset.

**NNET predicted curve using replicated dataset**


```r
#NNET
surv.predicted.nnet<-survfit(Surv(time, predicted.status)~ClusteredAge, data=nnet.test)
ggsurvplot(surv.predicted.nnet, conf.int = TRUE, censor= TRUE, cex.axis=3, cex.lab=3.0,  pval = TRUE, main="NNET predicted curve")
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-60-1.png)<!-- -->

**COX predicted curve using lung dataset**


```r
#COX
plot(tiempos1, c1.predicted.mean, ylim=c(0,1),col="indianred1",main="Predicted Survival curves",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(tiempos2,c2.predicted.mean,ylim=c(0,1),col="lightseagreen",main="Predicted Survival curves",xlab="Time",ylab="Survival probability", type="l", axes=FALSE)
legend("topright",c("Predicted (39,62] ","Predicted (62,82]"),lty=c(1,1),col=c("indianred1","lightseagreen"))
```

![](SurvivalAnalysis_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

Both curves are very similar to their reference one (the predicted nnet to the real using replicated data and the predicted cox to the real using lung dataset). 

The Cox predictions shows some peaks (as his respective real do) for the little amount of samples.


With the two predictions being so similar to real data for each dataset, I would say both have performed really well. ANN seems to be a little more accurate (as it's said the $P3-R2$ paper) but also needs more computing, time and data. 


