---
title: "Feature Selection"
author: "Naira María Chiclana García 44717497T"
date: "January 2019"
output: 
  html_document:
    self_contained: false
    keep_md: true
    toc: true
    toc_float: true
---



---

#Introduction to Feature Selection {#linkIntroduction}

Data can contain attributes that are highly correlated with each other. Many methods perform better if highly correlated attributes are removed. But, how can we choose the best subset to keep? 

**Feature Selection** is a optimization problem, and a component of Data Mining, specifically belonging to Data preprocessing (the first step, followed by modeling the data for make predictions and analysis of the results). Feature selection is done in order to **reduce the variables (features)** to work with (choosing an optimal or near-optimal subset with respect to an objetive function), as using all of them can lead to a bad performance and results.

Automatic feature selection methods can be used to build many models with different subsets of a dataset and identify those attributes that are and are not required to build an accurate model, and an adequate selection of features may improve accuracy and efficiency of classifier methods.

#### Strateggies:

1. **Search**: Examine all possible subsets and select the one that performs the best. The subsets will grow combinarotially (the selection of the optimal subset can't be guaranteed). 

    - **Optimun**: smallest possible search time

    - **Heuristic**:  optimize a problem by iteratively improving the solution based on a given heuristic function or a cost measure

    - **Randomized**: do not require the gradient of the problem to be optimized, and RS can hence be used on functions that are not continuous or differentiable



(2) **Evaluation**:


    - **Filter methods:** Selection based on general characteristics of the trainning set, it's independet of the classifier and it's fast and simple. Objetive function evaluates subsets by their information content. They have fast execution and generality, nevertheless, they have tendency to select large subsets.

    - **Wrapper methods:** Based on prediction accuracy of the chosen classifier, usually give better results but have a high computational cost. Objetive function evaluates subsets by their predictive accuracy (recognition rate) using statistical resampling or cross validation. The advantages of these methods are accuracy and ability to generalize while their disadvantages are slow execution and lack of generality.
    
-----

#Leukemia Dataset


From a collection of leukemia patient, 72 samples from bone marrow and peripheral blood. Contains measurements corresponding to $ALL$ (acute lymphoblast leukemia) 34 samples, and $AML$ (acute myeloid leukemia) 34 samples. *In acute disease the cells of bone marrow don't grow appropriately*.
We also have 7130 attributes defining each sample (as we can see below in the dimensions).



```r
ALL_data<-read.table("/Users/nairachiclana/Desktop/Computional Learning /Leo/P2-Feature selection/Leukemia files/Datos-Golub-Leuk/ALL.dat")
AML_data<-read.table("/Users/nairachiclana/Desktop/Computional Learning /Leo/P2-Feature selection/Leukemia files/Datos-Golub-Leuk/AML.dat")
```


```r
library(foreign)

ALL.AML_test<-read.arff("/Users/nairachiclana/Desktop/Computional Learning /Leo/P2-Feature selection/Leukemia files/ALL-AML_test.arff")
dim(ALL.AML_test)
```

```
## [1]   34 7130
```

```r
ALL.AML_train<-read.arff("/Users/nairachiclana/Desktop/Computional Learning /Leo/P2-Feature selection/Leukemia files/ALL-AML_test.arff")
dim(ALL.AML_train)
```

```
## [1]   34 7130
```

```r
cols<-ncol(ALL.AML_train)
myclass.train<-ALL.AML_train[,cols]
myclass.test<-ALL.AML_test[,ncol(ALL.AML_test)]
```

-----

#Lab 2-Part1

#### 1.1. Apply a filter and a wrapper method for selection of features from the provided data set.

### **Filter methods: ** 

As we have seen in [Introduction to Feature Selection](#linkIntroduction) this type of methods don't need a classifier, they carry out the feature selection process as a pre-processing step with no induction algorithm. And since it tends to select subsets with a high number of features (even all the features) a threshold will be  required to choose a subset (`subset.size`). We will use 10 by now (further in the analysis we will use  many more).

- **Chi-Squared: ** Makes use of *best.first.search*, a similar algorithm to *forward.search*, with the difference that this one chooses the best node of the already evaluated nodes and evaluated it.

- **Forest: ** Uses the qualifying *OneR* to find the weights of the attributes: for each attribute the algorithm create a rule based only on that attribute and then calculates the error ratio.




```r
library(FSelector)

aplicar.filtro <-function(data, filtro, subset.size) { 
  
  last.col <- colnames(data)[-1:-(ncol(data)-1)]
  formula<-as.formula(paste(last.col,".",sep="~"))
  
  #Chi-squared
  if(filtro=="chi") weigths<-chi.squared(formula, data)
  #Forest-importance
  else if(filtro=="forest") weigths<-random.forest.importance(formula, data)
  #10 best features
  best.weigths<-cutoff.k(weigths,subset.size)
  
  return(best.vbles=best.weigths)
}
```


```r
filter10.init<-Sys.time()
chi<-aplicar.filtro (ALL.AML_train, "chi", 10)
forest<-aplicar.filtro (ALL.AML_train, "forest", 10)
filter10.end<-Sys.time()
filter10.time<-filter10.end-filter10.init
```







```r
chi
```

```
##  [1] "attribute6855" "attribute758"  "attribute1685" "attribute1834"
##  [5] "attribute1882" "attribute2288" "attribute3252" "attribute1144"
##  [9] "attribute1902" "attribute2335"
```

```r
forest
```

```
##  [1] "attribute6855" "attribute3252" "attribute2288" "attribute1685"
##  [5] "attribute1882" "attribute2121" "attribute2335" "attribute1834"
##  [9] "attribute2361" "attribute1779"
```

```r
sum(chi %in% forest)
```

```
## [1] 7
```

```r
filter10.time
```

```
## Time difference of 51.53821 secs
```

7 of 10 attributes choosen by the different methods are the same (70% of coincidence).


```r
filter100.init<-Sys.time()
chi100<-aplicar.filtro (ALL.AML_train, "chi", 100)
forest100<-aplicar.filtro (ALL.AML_train, "forest", 100)
filter100.end<-Sys.time()
filter100.time<-filter100.end-filter100.init
```







```r
sum(chi100 %in% forest100)
```

```
## [1] 73
```

```r
filter100.time
```

```
## Time difference of 50.98928 secs
```

If we expand the search from 10 to 100, 73% of the choosen variables are the same, facing the 70% of coincidence with 10 variables, we could say that the bigger the subset, the more lightly similar the filters behave, but the coincidences tends to stay steady.

We can also see that the execution time is quite low, and  surprisingly, the time for 10 and 100 features are the same, illustrating the high efficiency  os this methods.




### **Wrapper method: **

**StepWise Forward Selection ** First, the best single feature is selected, and in each level n the algorithm computes which is the best combination of the best n features (the choosed ones with the remaining ones) basing on the *AUC* of the corresponding classifier. It will move forward for the *subset.size* (levels) and adding a feature (of the remaining ones) for each level, if the predefined number of features is already selected or result dont improve, the search will end.


 **Classifiers used: **

- **Decision Trees**: Logic construction diagrams are created to represent and categorize conditions which occur in a consecutive way. Parameters: *cp* are the branches (depth) and *minsplit* the minumun number of observations.

- **Neural Networks**: Basing o n artifical neurons connected among them, they form systems who learn and form themselves, re-establishing the neural weigths. Parameters: *size* is the number of neurons and *decay* to compensate the possible overfitting caused by too many neurons.

(The parameters have been previosly seelected by `caret`.)


```r
library(e1071)
library(rpart)
library(nnet)
library(plyr)
compute.models <-function(train.data, test.data, model, formula) {
  
  #Decision Trees
  if(model=="dt") {
   modelo<-rpart(formula,data=train.data,control=rpart.control(minsplit=3, cp=0.111));
   mod.estimado<-predict(modelo, newdata=test.data, type="matrix");
   mod.estimado<-as.matrix(mod.estimado)
  }
  
  #Neural network
  else if(model=="nnet") { 
    modelo<-nnet(formula,data=train.data, size=1, decay=0.1, maxit=3, trace=F);
    mod.estimado<-predict(modelo, newdata=test.data, type="raw");
  }
  
  #----------------Added to be used in CrossValidation function (Conclusions)-------
  #Logistic recresion
  else if(model=="glm") {
    modelo<-glm(formula, train.data, family=binomial("logit"));
    mod.estimado<-predict(modelo, newdata=test.data, type="response");
  }
  
  mod.estimado<-round(mod.estimado)
  ifelse(model =="dt",roc.test<-roc(as.numeric(test.data$myclass),as.numeric(mod.estimado[,dim(mod.estimado)[2]]),smooth=F,auc=T),roc.test<-roc(as.numeric(test.data$myclass), mod.estimado, smooth=F, auc=T)) 
  auc.test <-roc.test$auc
  return(list(auc=auc.test,roc=roc.test)) 
}
```



```r
library(FSelector)
library(caret)
library(pROC)

stepWise <- function(model, train.data, test.data, subset.size) {
  
  best.auc<-0
  auc.acumulado<-0
  lista.vbles<-c(colnames(train.data))
  lista.vbles <-lista.vbles[lista.vbles != "myclass"] 
  
  for(nivel in 1:length(lista.vbles)) { 
    best.auc<-0
    for(v in lista.vbles) { 
      nom_cols_comb<-v
      ifelse(nivel==1, att<-nom_cols_comb, att<-c(nom_cols_comb, lista.acumulada))
      ec <-as.simple.formula(att,"myclass")
      new.auc<-compute.models(train.data, test.data, model, ec)$auc[1];
      if(new.auc>best.auc) {
        best.auc<-new.auc
        best.modelo<-model
        best.vble.nivel<-nom_cols_comb
      }
    }
    lista.vbles<-lista.vbles[lista.vbles != best.vble.nivel]
    if(nivel!=1) auc.acumulado<-c(auc.acumulado, best.auc)
    ifelse(nivel==1, lista.acumulada<-best.vble.nivel,lista.acumulada<-c(lista.acumulada, best.vble.nivel))
    ifelse(nivel==1, auc.acum.ac<-best.auc, auc.acum.ac<-auc.acumulado[nivel-1])
    if(length(lista.acumulada)!=nivel || auc.acum.ac>best.auc[1] || length(lista.acumulada)==subset.size) break
  }
  best.formula<-as.simple.formula(lista.acumulada,"myclass")
  return(list(vbles=lista.acumulada, auc.acumulado=auc.acumulado, formula=best.formula))
}
```


```r
sw10.init<-Sys.time()
formula.sw.dt<-stepWise("dt",ALL.AML_train, ALL.AML_test, 10)
formula.sw.nnet<-stepWise("nnet",ALL.AML_train, ALL.AML_test, 10)
sw10.end<-Sys.time()
sw10.time<-sw10.end-sw10.init
```







```r
formula.sw.dt$vbles
```

```
##  [1] "attribute6855" "attribute1"    "attribute2"    "attribute3"   
##  [5] "attribute4"    "attribute5"    "attribute6"    "attribute7"   
##  [9] "attribute8"    "attribute9"
```

```r
formula.sw.nnet$vbles
```

```
## [1] "attribute1902" "attribute738"  "attribute142"  "attribute21"  
## [5] "attribute901"  "attribute653"  "attribute303"  "attribute5226"
## [9] "attribute828"
```

```r
sum(formula.sw.dt$vbles %in% formula.sw.nnet$vbles)
```

```
## [1] 0
```

```r
sw10.time
```

```
## Time difference of 17.79277 mins
```

There's not a single coincidence in the choosen attributes, and the execution time is **21.348 times bigger** than filters method's. 

---------------


#### 1.2. Select 100-200 features from the Leukemia data set using the selected filter and wrapper methods (for example correlation and SFS). We are going to use these subset of features for working with GAs.

We have already briefly described and seen how *Filter* and *Wrapper* methods perform. Let's see other different ones:


### Filter method: Correlation 


The correlation is a statistical tool used to measure the relationship between two or more variables. This measure can be used to analyze the correlation of a sample for one feature with the output (a high correlation value indicates that there is a relationship of the data with the output variable and thus we can select high correlated features). We can also have to correlated features that are also highly correlated between them. This might indicates redundancy of information, and in this case we can eliminate one of them.

Correlations uses method which evaluates subset of features on the basis of the hypothesis: *""Good feature subsets contain features highly correlated with the classification, yet uncorrelated to each other""*.

There are different correlation **methods** like

- **Pearson Correlation Coefficient (PCC) **: A number between [-1,1] indicates the extent to wich two variables are linearly related (Being 1 total positive linear correlation, 0 no linear correlation and -1 total negative linear correlarion).

The sample correlation coefficient is calculated by: 

$r=\frac{\sum(x_i-\bar{x})(y_i-\bar{y})}{\sqrt{\sum(x_i-\bar{x})^2} \sqrt{\sum(y_i-\bar{y})^2}}$

- **Spearman Correlation Coefficient **: "The Spearman correlation coefficient is defined as the Pearson correlation coefficient between the ranked variables.": Is equal to the Pearson correlation between the rank values of those two variables; while Pearson's correlation assesses linear relationships, Spearman's correlation assesses monotonic relationships (whether linear or not). If there are no repeated data values, a perfect Spearman correlation of +1 or −1 occurs when each of the variables is a perfect monotone function of the other.

Spearman rank formula:

$R_s=1-(\frac{\sigma \sum d^2}{n^3-n})$

- **Kendall rank correlation coefficient **: measures the ordinal association between two measured quantities using a *tau test* (a non-parametric hypothesis test for statistical dependence based on the tau coefficient). Follows the same rule of +1 and -1 that previous methods


Kendall $\tau$ coefficient:

- $\tau=\frac{number\ of \ concordant\ pairs - number\ of\ disconcordant\  pairs}{n(n-1)/2}$


**Functions **:

- `cor()` Computes the variance of features (all columns except the class).

- `findCorrelation()` searchs through a correlation matrix and returns a vector of integers corresponding to columns to remove to reduce pair-wise correlations. Parameter `cutoff` indicates the value for the pair-wise absolute correlation cutoff we will use 75%. 



```r
library(mlbench)
library(caret)
```

```
## Loading required package: lattice
```

```
## Loading required package: ggplot2
```

```r
features.correlation<-function(seed, data, subset.size, cut, method) { 
  
  set.seed(seed)
  #correlation matrix
  cols=ncol(data)
  correlationMatrix<-cor(data[,1:cols-1], method=method)
  #find attributes that are highly corrected 
  highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=cut)
  cutoff<-highlyCorrelated[1:subset.size]
  features.subset<-sapply(cutoff, function(i) colnames(data)[i])
  return(features.subset)
}
```


```r
correlation.subset.pearson<-features.correlation(7, ALL.AML_train, 100, 0.75, "pearson")
correlation.subset.spearman<-features.correlation(7, ALL.AML_train, 100, 0.75, "spearman")
correlation.subset.kendall<-features.correlation(7, ALL.AML_train, 100, 0.75, "kendall")
```






```r
correlation.subset.pearson
```

```
##   [1] "attribute21"  "attribute25"  "attribute34"  "attribute35" 
##   [5] "attribute46"  "attribute50"  "attribute70"  "attribute78" 
##   [9] "attribute94"  "attribute158" "attribute170" "attribute191"
##  [13] "attribute195" "attribute200" "attribute201" "attribute229"
##  [17] "attribute233" "attribute236" "attribute241" "attribute243"
##  [21] "attribute250" "attribute268" "attribute275" "attribute279"
##  [25] "attribute284" "attribute298" "attribute300" "attribute308"
##  [29] "attribute324" "attribute325" "attribute331" "attribute332"
##  [33] "attribute333" "attribute334" "attribute347" "attribute350"
##  [37] "attribute352" "attribute366" "attribute371" "attribute375"
##  [41] "attribute378" "attribute384" "attribute385" "attribute388"
##  [45] "attribute398" "attribute401" "attribute402" "attribute415"
##  [49] "attribute418" "attribute423" "attribute428" "attribute430"
##  [53] "attribute431" "attribute436" "attribute439" "attribute441"
##  [57] "attribute443" "attribute453" "attribute456" "attribute458"
##  [61] "attribute463" "attribute464" "attribute470" "attribute488"
##  [65] "attribute497" "attribute498" "attribute502" "attribute504"
##  [69] "attribute509" "attribute518" "attribute520" "attribute526"
##  [73] "attribute529" "attribute531" "attribute533" "attribute534"
##  [77] "attribute539" "attribute540" "attribute548" "attribute555"
##  [81] "attribute559" "attribute569" "attribute570" "attribute571"
##  [85] "attribute576" "attribute577" "attribute578" "attribute581"
##  [89] "attribute584" "attribute585" "attribute586" "attribute590"
##  [93] "attribute592" "attribute606" "attribute607" "attribute609"
##  [97] "attribute612" "attribute620" "attribute622" "attribute629"
```

```r
correlation.subset.spearman
```

```
##   [1] "attribute35"  "attribute40"  "attribute41"  "attribute94" 
##   [5] "attribute158" "attribute182" "attribute200" "attribute250"
##   [9] "attribute275" "attribute279" "attribute300" "attribute308"
##  [13] "attribute322" "attribute325" "attribute332" "attribute336"
##  [17] "attribute340" "attribute347" "attribute350" "attribute364"
##  [21] "attribute378" "attribute388" "attribute406" "attribute423"
##  [25] "attribute427" "attribute430" "attribute431" "attribute439"
##  [29] "attribute441" "attribute453" "attribute456" "attribute464"
##  [33] "attribute488" "attribute509" "attribute518" "attribute520"
##  [37] "attribute526" "attribute529" "attribute533" "attribute535"
##  [41] "attribute536" "attribute538" "attribute555" "attribute570"
##  [45] "attribute571" "attribute572" "attribute576" "attribute581"
##  [49] "attribute584" "attribute586" "attribute590" "attribute592"
##  [53] "attribute607" "attribute609" "attribute622" "attribute648"
##  [57] "attribute649" "attribute660" "attribute669" "attribute678"
##  [61] "attribute697" "attribute707" "attribute710" "attribute725"
##  [65] "attribute726" "attribute744" "attribute751" "attribute757"
##  [69] "attribute760" "attribute764" "attribute783" "attribute794"
##  [73] "attribute799" "attribute802" "attribute808" "attribute835"
##  [77] "attribute839" "attribute840" "attribute862" "attribute887"
##  [81] "attribute893" "attribute895" "attribute901" "attribute904"
##  [85] "attribute905" "attribute911" "attribute919" "attribute921"
##  [89] "attribute926" "attribute931" "attribute934" "attribute936"
##  [93] "attribute952" "attribute960" "attribute974" "attribute981"
##  [97] "attribute991" "attribute993" "attribute994" "attribute995"
```

```r
correlation.subset.kendall
```

```
##   [1] "attribute40"   "attribute1359" "attribute2026" "attribute3667"
##   [5] "attribute4652" "attribute4992" "attribute5058" "attribute5997"
##   [9] "attribute6224" "attribute6766" "attribute6803" "attribute6806"
##  [13] "attribute7005" "attribute7037" "attribute19"   "attribute5229"
##  [17] "attribute3836" NA              NA              NA             
##  [21] NA              NA              NA              NA             
##  [25] NA              NA              NA              NA             
##  [29] NA              NA              NA              NA             
##  [33] NA              NA              NA              NA             
##  [37] NA              NA              NA              NA             
##  [41] NA              NA              NA              NA             
##  [45] NA              NA              NA              NA             
##  [49] NA              NA              NA              NA             
##  [53] NA              NA              NA              NA             
##  [57] NA              NA              NA              NA             
##  [61] NA              NA              NA              NA             
##  [65] NA              NA              NA              NA             
##  [69] NA              NA              NA              NA             
##  [73] NA              NA              NA              NA             
##  [77] NA              NA              NA              NA             
##  [81] NA              NA              NA              NA             
##  [85] NA              NA              NA              NA             
##  [89] NA              NA              NA              NA             
##  [93] NA              NA              NA              NA             
##  [97] NA              NA              NA              NA
```


```r
sum(correlation.subset.pearson %in% correlation.subset.spearman)
```

```
## [1] 42
```

The coincidences are less than the half of the subset.


```r
library(dplyr)
train.correlation<-dplyr::select(ALL.AML_train, one_of(correlation.subset.pearson))
subset.CorrelationMatrix<-cor(train.correlation, method="spearman")
library(corrplot)
corrplot(subset.CorrelationMatrix, method="circle")
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-24-1.png)<!-- -->


```r
sum(correlation.subset.pearson %in% forest100)
```

```
## [1] 1
```

```r
sum(correlation.subset.pearson %in% chi100)
```

```
## [1] 1
```

```r
sum(correlation.subset.spearman %in% forest100)
```

```
## [1] 2
```

```r
sum(correlation.subset.spearman %in% chi100)
```

```
## [1] 3
```

We can barely see some rows and colums with a really good inverse relationship, but it's almost impossible to distinguish anything with so many variables. Let's see how these plots look for a smaller subset.

Between Pearson and Spearman, we will choose to visualize **Pearson Correlation** because the choosen features are more similar to the ones choosen by filter methods.



```r
correlation.subset10.spearman<-features.correlation(7, ALL.AML_train, 10, 0.75, "spearman")
train.correlation10<-dplyr::select(ALL.AML_train, one_of(correlation.subset10.spearman))
subset10.CorrelationMatrix<-cor(train.correlation10, method="spearman")
par(mfrow = c(1,2))
corrplot(subset10.CorrelationMatrix, method="circle")
corrplot(subset10.CorrelationMatrix, method="pie")
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

```r
par(mfrow = c(1,1))
corrplot(subset10.CorrelationMatrix, method="number")
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-26-2.png)<!-- -->
Obviously, the main diagonal is 1 because all features have positive total linear correlation with themselves. Also, it seems like there are more inverse than positive correlated variables. Almost all the pairs are around  $\simeq$  $\pm$ 0.3 correlation measures.


The **Kendall correlation**  only selected 17 attributes, let's see if this subset has better correlation than the others of 100.


```r
correlation.subset.kendall<-correlation.subset.kendall[1:17]
train.correlation.kendall<-dplyr::select(ALL.AML_train, one_of(correlation.subset.kendall))
subsetKendall.CorrelationMatrix<-cor(train.correlation.kendall, method="kendall")
corrplot(subsetKendall.CorrelationMatrix, method="circle")
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

The quality of the relations between the features don't have a remarkable improvement in front of Spearman, and pearson seems to be better, so **we are keeping Spearman for the rest of the analysis**. Another thing that our attention is caught by, is that there arent any pair of features with total negative linear correlarion or less than -0.4, while in Spearman there are some of them.


### Wrapper method: Backwards Feature Selection (SBS)  

As we have said and seen before, wrapper methods needs a classifier whose prediction accuracy is the basis.

Backwards feature elimination is the same as recursive feature elimination (**RFE**): First, the critetion function is computed for all the features in the dataset (the algorithm fits the model to all predictors), and through the irerations the features are deleted and criterion function is computed for all subsets with *n-i* features, discarting the wrost one, until only the predefined number of features are left.

As a classifier, we will use a control object with 10 folds cross-validation (made with function `rfeControl()`). The function for model fitting is one predefined in the package `mlbench`: `rfFuncs`.

With `rfe()` (the function whose execute the RFE) the predictors (features) are ranked and the less important ones are sequentially eliminated. How the `size` parameter affects: t each iteration of feature selection, the $S_{i}$ top ranked predictors are retained, the model is refit and performance is assessed. The value of $S_{i}$ with the best performance is determined and the top $S_{i}$ predictors are used to fit the final model. 


```r
library(caret)
library(mlbench)
```


```r
features.backwards<-function(data, subset.size) {
  rfe_controller <- rfeControl(functions=rfFuncs, method="cv", repeats=10)
  rfe.subset<-rfe(data[,1:cols-1], data[,cols], rfeControl=rfe_controller, size=100)
  features.subset<-predictors(rfe.subset)
  return(features.subset)
}
```


```r
sbs.init<-Sys.time()
sbs.subset<-features.backwards(ALL.AML_train, 100)
sbs.end<-Sys.time()
sbs100.time<-sbs.end-sbs.init
```







```r
sbs100.time
```

```
## Time difference of 22.02901 secs
```

```r
sbs.subset
```

```
##   [1] "attribute6855" "attribute1834" "attribute6041" "attribute3252"
##   [5] "attribute1685" "attribute758"  "attribute1882" "attribute2288"
##   [9] "attribute1144" "attribute2642" "attribute1829" "attribute1902"
##  [13] "attribute4847" "attribute760"  "attribute3469" "attribute2141"
##  [17] "attribute4894" "attribute4204" "attribute2361" "attribute4644"
##  [21] "attribute2335" "attribute2128" "attribute2354" "attribute6225"
##  [25] "attribute5171" "attribute4366" "attribute2121" "attribute6376"
##  [29] "attribute4680" "attribute1779" "attribute6725" "attribute5501"
##  [33] "attribute4229" "attribute1725" "attribute5688" "attribute7119"
##  [37] "attribute922"  "attribute4463" "attribute4116" "attribute2363"
##  [41] "attribute6919" "attribute2714" "attribute2402" "attribute3897"
##  [45] "attribute804"  "attribute5280" "attribute2422" "attribute6522"
##  [49] "attribute5300" "attribute3004" "attribute4377" "attribute4318"
##  [53] "attribute4196" "attribute1159" "attribute4925" "attribute6347"
##  [57] "attribute4373" "attribute6185" "attribute671"  "attribute6257"
##  [61] "attribute6005" "attribute3804" "attribute5833" "attribute6761"
##  [65] "attribute354"  "attribute235"  "attribute3482" "attribute2497"
##  [69] "attribute1369" "attribute2441" "attribute6573" "attribute1909"
##  [73] "attribute5955" "attribute1953" "attribute3475" "attribute4763"
##  [77] "attribute4082" "attribute1157" "attribute1598" "attribute6250"
##  [81] "attribute1400" "attribute4480" "attribute4328" "attribute4008"
##  [85] "attribute88"   "attribute3904" "attribute6281" "attribute3775"
##  [89] "attribute5335" "attribute3104" "attribute885"  "attribute332" 
##  [93] "attribute3422" "attribute6542" "attribute2014" "attribute5377"
##  [97] "attribute1674" "attribute3237" "attribute2833" "attribute4582"
```

Surprisingly, we were expecting a high execution time since it's a wrapper method, but this one is even faster than filter methods.

------------------


#Introduction to Genetic Algorithms (GAs)

Genetic Algorithms are search and optimization techniques based on Darwin’s Principle of Natural Selection.

Inspired by the biological mechanism of natural selection and reproduction, this algorithm has a global optimization technique for searching very large spaces. It uses a probabilistic search using a population of encoded structure. Uses objective function. The possible solutions are represented as chromosomes (finite sequences of binary codification).

**Search and genetic algorithem procces:**

(1) **Inizialization**: Determine the number of chromosomes, generation, and mutation rate and crossover rate value. Generate chromosome-chromosome number of the population, and the initialization value of the genes chromosome-chromosome with a random value.

(2) **Evaluation** of fitness value of chromosomes (offspring population) by calculating objective function.

(3) **Selection** of chromosomes. It will save copies of the solutions with higher fitness values (explained in [Fitness function](#linkfitness) ).

(4) **Recombination (Crossover)**  of parts of two or more parental solution to create new and possible better solutions.

(5) **Mutation** locally but randomly in solutions.

(6) **Replacement** of the old chromosomes (parents) by the new offsprings (the ones created, recombinated and mutated). 

(7) Repeat from Evaluation till replacement until the therminating condition is met, we prpbably end with the best chromosomes.




We will to perform **binary encode** to the solution, where chromosomes are strings of 1s and 0s and each position in the chromosome represents a particular characteristic of the problem.

-----------------

#Lab 2-Part2

#### 2.1. In order to implement a Genetic algorithm, a prediction function is needed in order to evaluate the fitness function. Use a simple predictor function like LDA.


**Linear discriminant analysis (LDA)** :  

The predictor function LDA uses linear combinations of predictors to predict the class of a given observation. 

LDA assumes that predictors are normally distributed (Gaussian distribution) and that the different classes have class-specific means and equal variance/covariance. 


```r
library(MASS)
library(tidyverse)
```

How it works: 

The algorithm starts by finding directions that maximize the separation between classes (`lda()`) , then use these directions to predict the class of individuals (`predict()`). These directions, called linear discriminants, are a linear combinations of predictor variables.


```r
# Fit the model
model<-lda(myclass~., data=ALL.AML_train)
# Make predictions
predictions <- model %>% predict(ALL.AML_test)
# Model accuracy
mean(predictions$class==ALL.AML_test$myclass)
```

```
## [1] 0.8529412
```


------------

###Fitness evaluation {#linkfitness}

Performs by a **Fitness function** who quantifies the optimality of a solution (chromosome) so that that particular solution may be ranked against all the other solutions. A **fitness value** is assigned to each solution depending on **how close it actually is to solving the problem**. 
Ideal fitness function correlates closely to goal and is quickly computable.

The used formula will be: 
$fitnes=-C_{1}*error-C_{2}*sum(ones)*$ 
Favoring solutions with lower prediction error and penalising a higher number of selecter attributes.


We will make one function with several parameters that will be used to best combination  of:

- Propbability of 1 in the parent chromosome.

- Factor value (used in penalty).



```r
library(MASS)
library(FSelector)
library(mclust)
```


```r
fitness<-function(train, test, prob, penalty, factor) {
  cols<-ncol(train)
  indices<-as.vector(rbinom(cols,1, prob))
  indices[cols]<-0
  result=-1
  if(sum(indices)>1) {
    model<-lda(as.simple.formula(colnames(train)[indices==1], "myclass"),train)
    predictions<-model %>% predict(test)
    result=-classError(predictions$class, test$myclass)$errorRate
    if(penalty==T) result=result-factor*sum(indices)/500
  }
  return(result)
}
```

The result of this function is the error, the ideal value is the nearest possible to 0.


```r
#Witouth penalty
init005<-Sys.time()
fitness(ALL.AML_train, ALL.AML_test, 0.05, F, 0.8)
```

```
## [1] -0.1176471
```

```r
end005<-Sys.time()
time005<-end005-init005
time005
```

```
## Time difference of 0.1052461 secs
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.1, F, 0.8)
```

```
## [1] -0.1176471
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.3, F, 0.8)
```

```
## [1] -0.1470588
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.5, F, 0.8)
```

```
## [1] -0.1176471
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.7, F, 0.8)
```

```
## [1] -0.1470588
```

```r
init09<-Sys.time()
fitness(ALL.AML_train, ALL.AML_test, 0.9, F, 0.8)
```

```
## [1] -0.1176471
```

```r
end09<-Sys.time()
time09<-end09-init09
time09
```

```
## Time difference of 20.91022 secs
```

By general terms, it seems like the more the variables, the less the error when no penalty is used. The computing time increases a lot when addig more features. Probably if we use penalty with will increase with the increase of variables. Lets try the best ones to choose a suitable probability for the input chromosome.



```r
#With penalty
fitness(ALL.AML_train, ALL.AML_test, 0.1, T, 0.8)
```

```
## [1] -1.332282
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.3, T, 0.8)
```

```
## [1] -3.686259
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.5, T, 0.8)
```

```
## [1] -5.908659
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.9, T, 0.8)
```

```
## [1] -10.46866
```

Effectively, it is worse when there is more features, we will chose 0.3 because our subset is very large, and it's better tho choose less but significant variables  (lower number of features and higher penalty).


```r
#Changing factor
fitness(ALL.AML_train, ALL.AML_test, 0.3, T, 0.4)
```

```
## [1] -1.842259
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.3, T, 1)
```

```
## [1] -4.369647
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.6, T, 0.4)
```

```
## [1] -3.545647
```

```r
fitness(ALL.AML_train, ALL.AML_test, 0.6, T, 1)
```

```
## [1] -8.701059
```

Obviously, the error is lower with a lower penalty, but the results are less reliable to be the optimun result.

A medium rare could be:


```r
fitness(ALL.AML_train, ALL.AML_test, 0.6, T, 0.7)
```

```
## Warning in lda.default(x, grouping, ...): variables are collinear
```

```
## [1] -6.164247
```

Once we have find the more suitable parameters of probability of 1 in chromosome ($0.6$)  and  penalty ($0.7$) we make two functions with these ones fixed that will be called by `ga()` algorithm with the choosen subsets. One with  penalty `fitness.with.penalty.ga()` and another without it `fitness.ga()`.


```r
fitness.ga<-function(train, test) {
  result=fitness(train,test, 0.6, F, 0.7)
  return(result)
}

fitness.with.penalty.ga<-function(train, test) {
  result=fitness(train,test, 0.6, T, 0.7)
  return(result)
}
```


----


#### 2.2.  Use a filter and/or a wrapper method (Correlation/SFS) to reduce the number of input variables from 7129 to a number that permits you to run the GA algorithm.

`ga()` function performs a maximization of our [Fitness function](#linkfitness) using genetic algorithms. Local search using general- purpose optimisation algorithms can be applied stochastically to exploit interesting regions. 

### GA algorithm {#linkGA}


```r
library(GA)
```


```r
#All parameters
init.all<-Sys.time()
GA.all<-ga("binary", fitness= function(x) fitness.ga(ALL.AML_train, ALL.AML_test),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter = 180)
end.all<-Sys.time()
time.all<-end.all-init.all
```






If we use all the parameters in the dataset, the execution time is 


```r
time.all
```

```
## Time difference of 8.163058 hours
```


It is clear that we can use all the parameters (7130) because it's extremely inefficient.


Let's use some of the filter and wrapper methods shown above to make different sizes subsets till we find a suitable size for an efficient execution.


**1000 features** with Filter methods (the fastest ones):


```r
#Chi subsets
library(dplyr)
subset1000.filter.chi<-aplicar.filtro (ALL.AML_train, "chi", 1000)
train1000chi<-dplyr::select(ALL.AML_train, one_of(subset1000.filter.chi))
test1000chi<-dplyr::select(ALL.AML_test, one_of(subset1000.filter.chi))
#add last colum if it has been removed by filter
if(subset1000.filter.chi[1000] != "myclass")  {
  train1000chi$myclass<-ALL.AML_train$myclass
  test1000chi$myclass<-ALL.AML_test$myclass
}

  
#chi ga
init.1000.chi<-Sys.time()
GA.1000.chi<-ga("binary", fitness= function(x) fitness.ga(train1000chi, test1000chi),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter = 180)
end.1000.chi<-Sys.time()
time.1000.chi<-end.1000.chi-init.1000.chi

#Forest subsets
subset1000.filter.forest<-aplicar.filtro (ALL.AML_train, "forest", 1000)
train1000forest<-dplyr::select(ALL.AML_train, one_of(subset1000.filter.forest))
test1000forest<-dplyr::select(ALL.AML_test, one_of(subset1000.filter.forest))
#add last colum si la ha quitado con filtro
if(subset1000.filter.forest[1000] != "myclass")  {
  train1000forest$myclass<-ALL.AML_train$myclass
  test1000forest$myclass<-ALL.AML_test$myclass
}
#forest ga
init.1000.forest<-Sys.time()
GA.1000.forest<-ga("binary", fitness= function(x) fitness.ga(train1000forest, test1000forest),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter = 180)
end.1000.forest<-Sys.time()
time.1000.forest<-end.1000.forest-init.1000.forest
```



```r
time.1000.chi
```

```
## Time difference of 10.16857 mins
```

```r
time.1000.forest
```

```
## Time difference of 10.37428 mins
```

Both have more or less the same time, we will choose *Chi-Squared* which barely beat *Random Forest*.

**500 features** with Chi-Squared Filter: 


```r
#Forest subsets
subset500.filter.chi<-aplicar.filtro (ALL.AML_train, "chi", 500)
train500chi<-dplyr::select(ALL.AML_train, one_of(subset500.filter.chi))
test500chi<-dplyr::select(ALL.AML_test, one_of(subset500.filter.chi))
#add last colum si la ha quitado con filtro
if(subset500.filter.chi[500] != "myclass")  {
  train500chi$myclass<-ALL.AML_train$myclass
  test500chi$myclass<-ALL.AML_test$myclass
}
#forest ga
init.500.chi<-Sys.time()
GA.500.chi<-ga("binary", fitness= function(x) fitness.ga(train500chi, test500chi),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter = 180)
end.500.chi<-Sys.time()
time.500.chi<-end.500.chi-init.500.chi
```



```r
time.500.chi
```

```
## Time difference of 3.536188 mins
```

Reducing the subset to the half, the time reduces more than the half. 3 mins appears to be a reasonable execution time.
We will also try with the same size of subset but the variables choosen by the **correlation** method and see if it goes faster.



```r
subset500.correlation<-features.correlation(7, ALL.AML_train, 500, 0.75, "spearman")
train500correlation<-dplyr::select(ALL.AML_train, one_of(subset500.correlation))
test500correlation<-dplyr::select(ALL.AML_test, one_of(subset500.correlation))
#add last colum si la ha quitado con filtro
if(subset500.correlation[500] != "myclass")  {
  train500correlation$myclass<-ALL.AML_train$myclass
  test500correlation$myclass<-ALL.AML_test$myclass
}
#forest ga
init.500.correlation<-Sys.time()
GA.500.correlation<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter = 180)
end.500.correlation<-Sys.time()
time.500.correlation<-end.500.correlation-init.500.correlation
```



```r
time.500.correlation
```

```
## Time difference of 3.718099 mins
```

Chi filter keeps in the first place regarding to the execution time.

Let's also see the subset choosen by **SBS**:


```r
subset500.sbs<-features.backwards(ALL.AML_train, 500)
train.subset.sbs<-dplyr::select(ALL.AML_train, one_of(subset500.sbs))
test.subset.sbs<-dplyr::select(ALL.AML_test, one_of(subset500.sbs))
#add last colum si la ha quitado con filtro
if(subset.sbs[length(subset.sbs)] != "myclass")  {
  train.subset.sbs$myclass<-ALL.AML_train$myclass
  test.subset.sbs$myclass<-ALL.AML_test$myclass
}
#forest ga
init.sbs<-Sys.time()
GA.500.sbs<-ga("binary", fitness= function(x) fitness.ga(train.subset.sbs, test.subset.sbs),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter = 180)
end.sbs<-Sys.time()
time.500.sbs<-end.sbs-init.sbs
```





```r
time.500.sbs
```

```
## Time difference of 48.70931 secs
```

This time is BY FAR the fastest one optimizing the fitness function. We will later see if the *AUC* resulst are as good.

We will visualize all the plots with the same interval in Y axis to see the differences of fitness values clearly.




```r
#all features
#summary(GA.all)
#plot(GA.all, ylim=c(0, 0.3))

#chi 1000 subset
summary(GA.1000.chi)
```

```
## Loading required package: GA
```

```
## Loading required package: foreach
```

```
## 
## Attaching package: 'foreach'
```

```
## The following objects are masked from 'package:purrr':
## 
##     accumulate, when
```

```
## Loading required package: iterators
```

```
## Package 'GA' version 3.2
## Type 'citation("GA")' for citing this R package in publications.
```

```
## 
## Attaching package: 'GA'
```

```
## The following object is masked from 'package:utils':
## 
##     de
```

```
## ── Genetic Algorithm ─────────────────── 
## 
## GA settings: 
## Type                  =  binary 
## Population size       =  30 
## Number of generations =  180 
## Elitism               =  2 
## Crossover probability =  0.8 
## Mutation probability  =  0.1 
## 
## GA results: 
## Iterations             = 180 
## Fitness function value = -0.02941176 
## Solutions = 
##      x1 x2 x3 x4 x5 x6 x7 x8 x9 x10  ...  x499 x500
## [1,]  0  1  0  0  0  0  1  0  0   1          1    1
## [2,]  0  1  0  0  0  0  1  0  0   1          1    1
```

```r
plot(GA.1000.chi, ylim=c(-0.08, 0))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-58-1.png)<!-- -->

```r
#chi 500 subset
summary(GA.500.chi)
```

```
## ── Genetic Algorithm ─────────────────── 
## 
## GA settings: 
## Type                  =  binary 
## Population size       =  30 
## Number of generations =  180 
## Elitism               =  2 
## Crossover probability =  0.8 
## Mutation probability  =  0.1 
## 
## GA results: 
## Iterations             = 180 
## Fitness function value = -0.02941176 
## Solutions = 
##       x1 x2 x3 x4 x5 x6 x7 x8 x9 x10  ...  x499 x500
## [1,]   0  1  0  1  1  0  0  1  1   1          0    0
## [2,]   0  1  0  1  1  0  0  1  1   1          1    0
## [3,]   0  1  0  1  1  0  0  1  1   1          0    1
## [4,]   0  1  0  1  1  0  0  1  1   0          0    1
## [5,]   0  1  0  1  1  0  0  1  1   0          0    1
## [6,]   0  1  0  1  1  0  0  1  1   0          0    1
## [7,]   0  1  0  1  1  0  0  1  1   0          0    1
## [8,]   0  1  0  1  1  0  0  1  1   0          0    0
## [9,]   0  1  0  1  1  0  0  1  1   1          0    0
## [10,]  0  1  0  1  1  0  0  1  1   1          0    1
##  ...                                                
## [16,]  0  1  0  1  1  0  0  1  0   0          0    1
## [17,]  0  1  0  1  1  0  0  1  1   0          0    1
```

```r
plot(GA.500.chi, ylim=c(-0.08, 0))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-58-2.png)<!-- -->

```r
#correlation subset
summary(GA.500.correlation)
```

```
## ── Genetic Algorithm ─────────────────── 
## 
## GA settings: 
## Type                  =  binary 
## Population size       =  30 
## Number of generations =  180 
## Elitism               =  2 
## Crossover probability =  0.8 
## Mutation probability  =  0.1 
## 
## GA results: 
## Iterations             = 180 
## Fitness function value = -0.02941176 
## Solutions = 
##      x1 x2 x3 x4 x5 x6 x7 x8 x9 x10  ...  x499 x500
## [1,]  0  1  0  0  1  0  0  0  1   1          1    1
## [2,]  0  1  0  0  1  0  0  0  1   1          1    1
## [3,]  0  1  0  0  1  0  0  0  1   1          1    1
```

```r
plot(GA.500.correlation, ylim=c(-0.08, 0))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-58-3.png)<!-- -->

```r
#SBS subset
summary(GA.500.sbs)
```

```
## ── Genetic Algorithm ─────────────────── 
## 
## GA settings: 
## Type                  =  binary 
## Population size       =  30 
## Number of generations =  180 
## Elitism               =  2 
## Crossover probability =  0.8 
## Mutation probability  =  0.1 
## 
## GA results: 
## Iterations             = 180 
## Fitness function value = 0 
## Solutions = 
##       x1 x2 x3 x4 x5 x6 x7 x8 x9 x10  ...  x499 x500
##  [1,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [2,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [3,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [4,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [5,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [6,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [7,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [8,]  1  0  0  1  1  0  0  1  1   1          1    1
##  [9,]  1  0  0  1  1  0  0  1  1   1          1    1
## [10,]  1  0  0  1  1  0  0  1  1   1          1    1
## [11,]  1  0  0  1  1  0  0  1  1   1          1    1
## [12,]  1  0  0  1  1  0  0  1  1   1          1    1
```

```r
plot(GA.500.sbs,ylim=c(-0.08, 0))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-58-4.png)<!-- -->

For all the subsets we have used the fixed parameters:

- Probabilty of crossover between chromosomes=0.8

- Probability of mutation in a parent chromosome=0.1

- Number of bits used in the binary encoded=500

- Population size=30

- Number of generations=180

And the results of the fitness value function and the execution time of ga have been:

| Subset used  | Execution time | Mean fitness value around |
|-------------|----------------|
|All parameters  | 8 hours  |    |
|Filter  Chi 1000 |  10.16 mins   |  -0.075  |
|Filter  Chi 500 |   3.54 mins  |  -0.03 |
|Filter Correlation 500 |  3.72 mins |   -0.07  |
|Wrapper SBS |  48.71 secs  |  -0.01  |

Since we're visualizing the error, the best value is the one nearest to 0.

The best result fitness value is the same for all subsets. Nevertheless, the one with the best mean is SBS, with all the features around -0.01 value, and the execution time is by far the shortest one. This result seems strange, for this reason, we're keeping the second best one too: Chi500, with the values around -0.03.

Chi with 100 features and correlation with 500 are very similar, but 1000 features chi has lower peaks. **Here we have the proof: more variables not only increase the execution time but run a worse performance.** of course, we're going to discard the 1000 variables subset.

As we have said, correlation is such a good result that is hard to believe, consequenly we're going to check the AUC to see if it's as good as it seems. (Function `crossValidation()` in [Conclusions](#linkConclusions) )



```r
formula.sbs500<-as.simple.formula(paste(subset500.sbs, collapse="+"), "myclass")
library(pROC)
set.seed(7)
all.data.sbs<-rbind(train.subset.sbs, test.subset.sbs)
aucs<-crossValidation(all.data.sbs, 10, formula.sbs500);
```

```r
aucs
```

```
## $glm.mean
## [1] 0.9416667
## 
## $dt.mean
## [1] 1
## 
## $nnet.mean
## [1] 0.6041667
```


I's normal that Neural Networks don't give a very good result because this classifier need a big amount of data to train the neurons, and we have keep only a "small" subset of 500. Anyways, the result is good, and with Logistic regression and Decision trees the AUC is almost perfet. With no doubt, **we're choosing and keeping Selective Backwards Selection subset**.


#### 2.3. Play with the parameters of the GA function to optimize operation of the genetic algorithm. Use a penalty term for reducing the number of features, and analyze that is indeed working.

### GA parameters and penalty {#link23}

Since the result of SBS is almost unbeatable, for seeing the change with parameters we're going to use another improvable subset: Correlation

Some of the avaliable parameters for this complex fuctions and that we will work with:

- `type` the type of genetic algorithm to be run depending on the nature of decision variables. Here we will alwaus use "binary" for binary representations of decision variables.

- `fitness` the fitness function, any allowable R function which takes as input an individual string representing a potential solution, and returns a numerical value describ- ing its “fitness”. We will use [our fitness function](#linkfitness)

- `nBits` a value specifying the number of bits to be used in binary encoded optimizations. We will keep using the subset of 500 variables, for this reason this parameter won't neither change.

- `maxiter` the maximum number of iterations to run before the GA search is halted (number of generations).


- `elistism` the number of best fitness individuals to survive at each generation. By default the top 5% individuals will survive at each iteration.

- `popSize` the population size.

- `pcrossover` the probability of crossover between pairs of chromosomes. Typically this is a large value and by default is set to 0.8.

- `pmutation` the probability of mutation in a parent chromosome. Usually mutation occurs with a small probability, and by default is set to 0.1.

- `run` the number of consecutive generations without any improvement in the best fitness value before the GA is stopped.




It is self-evident that the line refering to the best value don't need that much **iterations** to get steady value. With all the variables have a slop in 70, in forest and SBS doesn't change since the begining, and in correlation the slop is more or less at iteration 10. So, more than 100 iterations are unnecessary when using all variables (never), while for a subset of 500 is enough with 50.

Anyway, let's confirm that the result dosn't get worse when reducing the generations from 180 to 50.



```r
GA.correlation.50it<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter=50)
```







```r
plot(GA.correlation.50it, ylim=c(-0.09, -0.02))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-65-1.png)<!-- -->

```r
time.correlation.50it
```

```
## Time difference of 1.058009 mins
```

The fitness value has has remained in the same value and the time is three times less. From now on, we're using 50 iterations instead of 180, the value hasn't improved but the execution has.

Regarding to **elistism** value, we expect that if we keep more of the best individuals for sure, the result well be better, let's ascertain this assumption, for an exaggerated example, we will try to kepp the 20%. (It was 5% by default)


```r
GA.correlation.bigger.elitism<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter=50, elitism=0.2)
```




```r
plot(GA.correlation.bigger.elitism,,  ylim=c(-0.09, -0.02))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-67-1.png)<!-- -->

The result is s little bit worse, barely new individuals (possible better) have been picked and the error is bigger.

In relation to **run** parameter, we expect that if we reduce the threshold of generations without any improvement allowed to continue, the result will improve. Let's check that reasoning:


```r
GA.correlation.lower.run<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter=80, elitism=0.2, run=2)
```





```r
plot(GA.correlation.lower.run)
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-69-1.png)<!-- -->
The result is worse, probably because there were no enough iterations, let`s try with `run=10`:


```r
GA.correlation.run10<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.8, pmutation=0.1, nBits = 500, monitor=FALSE, popSize = 30, maxiter=80, elitism=0.2, run=10)
```




```r
plot(GA.correlation.run10, ylim=c(-0.09, -0.02))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-71-1.png)<!-- -->

If it's to low, the lack of iterations remain the value really slow, and a reasonable one (10 of 80) doesn't neither have a remarkable difference, in any case, it's curious how in such a low amount of iterations (4) the results reach the same value that all the previous ones.

The parameters **pcrossover** and **pmutation** provides us nothing but diversity and randomness which could be either good or bad, we can be sure of what is going to happen but let's see:

Increase the mutation and reduce the crossover:


```r
GA.mutation.bigger<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.2, pmutation=0.6, nBits = 500, monitor=FALSE, popSize = 30, maxiter=80, elitism=0.2)
```

Decrease the mutation and increase the crossover:


```r
GA.crossover.bigger<-ga("binary", fitness= function(x) fitness.ga(train500correlation, test500correlation),  pcrossover=0.95, pmutation=0.01, nBits = 500, monitor=FALSE, popSize = 30, maxiter=80, elitism=0.2)
```




```r
plot(GA.crossover.bigger, ylim=c(-0.09, -0.02))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-74-1.png)<!-- -->

```r
plot(GA.mutation.bigger,  ylim=c(-0.09, -0.02))
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-74-2.png)<!-- -->

In spite of the changes have been really drastic, the value has neither changed significantly, being a bigger crossover a better option than a bigger % of mutations.

Best parameters performance: **Reduce iterations** to mantain fitness value but reduce execution time and  **increase the crossovers**.



What happens when we also use a **penalty** for keeping more variables?


```r
GA.with.penalty<-ga("binary", fitness= function(x) fitness.with.penalty.ga(train500correlation, test500correlation),  pcrossover=0.5, pmutation=0.3, nBits = 500, monitor=FALSE, popSize = 30, maxiter=50)
```


```r
summ.penalty<-summary(GA.with.penalty)
summ.penalty$fitness
```

```
## [1] -0.4200118
```

```r
plot(GA.with.penalty)
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-76-1.png)<!-- -->

The result has gotten much worse, hat was to be expected since it's the first time we're decreasing the value even more for the features picked. Further on, we will see if the AUC value is better or worse.




-----


#### 2.4.Analyze different executions of the algorithms to see if the same variables are obtained in each run, or at least which are the most frequent ones. {#link24}

We are going to use the `summary` parameter `$solution`, a matrix of solutions strings with as many rows as the number of solutions founds and as many columns as the number of decision variables.

For each one of the results given by the in the differents executions of `ga()` algorithm (in section [Exploration of parameters in GA function](#link23)), we extract the variables of the subset used in ga (*subset500.correlation*) that are 1 in the solution matrix, that is, the ones choosen for the final population solution. 


```r
features.frecuency<-function(ga.results.list) {
  list.of.solution.features<-list()
  features.in.solution<-list()
  for(ga.result in 1:length(ga.results.list)) {
    summ<-summary(ga.results.list[[ga.result]])
    solution.matrix<-summ$solution
    for(col in 1:ncol(solution.matrix)) {
      if(solution.matrix[1,col]==1) {
        features.in.solution<-c(features.in.solution, subset500.correlation[col])
      }
    }
    features.in.solution<-unlist(features.in.solution)
    list.of.solution.features[[ga.result]]<- features.in.solution
    
  }
  return(list.of.solution.features)
}
ga.results<-c(GA.with.penalty,GA.correlation.50it, GA.correlation.bigger.elitism, GA.correlation.lower.run,GA.correlation.run10, GA.mutation.bigger, GA.crossover.bigger)
list.appearances<-features.frecuency(ga.results)
ranking<-table(unlist(list.appearances))
sort(ranking, decreasing=T)[1:50]
```

```
## 
## attribute1621  attribute757 attribute2235 attribute2335 attribute2871 
##            28            28            27            27            27 
## attribute3246 attribute1835  attribute275 attribute2809 attribute3718 
##            27            26            26            26            26 
## attribute3763 attribute3849  attribute726 attribute1084 attribute1812 
##            26            26            26            25            25 
## attribute1882 attribute2002 attribute2599 attribute3040 attribute3341 
##            25            25            25            25            25 
## attribute3870 attribute1000 attribute1224 attribute2611 attribute2855 
##            25            24            24            24            24 
## attribute3778  attribute533 attribute1257 attribute1551 attribute1745 
##            24            24            23            23            23 
## attribute2322 attribute2508 attribute2743 attribute3882  attribute406 
##            23            23            23            23            23 
##  attribute584 attribute1110 attribute1998 attribute2462 attribute2565 
##            23            22            22            22            22 
## attribute2667 attribute2708 attribute2742 attribute2848 attribute3242 
##            22            22            22            22            22 
##   attribute35 attribute3707 attribute3738  attribute744  attribute760 
##            22            22            22            22            22
```
 
*There are more counts that the number of solutions because inside the same population, the chromosome with the best fitness don't have to be unique.*

------

### Conclusions and AUCs {#linkConclusions}

**Execution time** of all the methods used and the size of subsets choosen:

| Method | Subset size | Execution time | 
|--------|--------|--------------------|
|Filter  Chi Squared |  10 |  49 secs   |  
|Filter Random Forest | 10 | 51 secs |   
|Handmade Wrapper |  10 | 9.4 mins | 
|Filter Correlation |100 |  36 secs|
|Filter Chi Squared  |100 |  50 secs| 
|Filter  Random Forest | 100 |   50 secs |  
|SBS Wrapper |100 |  22 secs| 
|SBS Wrapper |500 | 49 secs |
|Filter Chi Squared | 500 |   3.5 mins |
|Filter Correlation | 500 |   3.7 mins |
|Filter Chi Squared | 1000 |   10.2 mins |
|Filter Random Forest | 1000 |  10.37 mins |




In [Introduction to Feature Selection](#linkIntroduction) we said that Filter methods are much more faster than Wrapper methods. With the handmade SFS algorihtm, the difference of time is very notable, more than 10 times slower than Filters, but, surprisingly, the used **SBS Wrapper** method from `mlbench` library, is, with a big diffecence, **the fastest** of all the ones used in all the study. The **slowest** is the **handmade SFS Wrapper** (as expected). Filters Chi and Forest don't have a notable difference, being Chi slightly fastest for larger subsets. The Filter method which used the less execution time was Correlation.

For the [Fitness Function](#linkfitness) used in GA algorithm, following several tests we concluded that the most suitable parameters were $60\%$ of probability of 1 in the input chromosome and a penalty of $0.7*sum(indices)/500$.

When testing the **most appropiate number of features** for [GA Algorithm](#linkGA), we saw that using all the parameters, not only the execution is unviable, but the results were worse that with a smaller subset, **proving that Feature Selection is essential with working with large datasets**. We choosed a final subset size of 500 features, since the execution time was not bigger than 3 mins, and the GA algorithm would 	hypothetically improve the results. 

The methods used for choosing the input subsets were: Chi-Squared, SBS and Correlation. All the results had a best value reached of -0.029 and the mean around -0.07 (which is pretty good), but to our surprise, the **subset provided by SBS** (the fastest one) was, with a major difference, **the best one**, reaching a best value of 0 (perfect) and with the mean around -0.01. The result was so good that was hard to believe, so we test the AUC with the classifiers *Linear Regresion, Decision Trees and Neural Networks* to verify if it really was that good and it seemed. Unvelievably the AUC confirmed the theory, giving an AUC of $1$ and $0.94$.

For [Parameters exploration](#link23) we used a worse subset (the one given by correlation) since SBS was unbeatable. We find that, for this particular dataset: That much iterations (180) were uneeded since the values got steady at more or less 20, and the computational cost was decreased. Increase the number of offspring keeped (*elitism*) deteriorate the results, which seems to indicate that "evolution" is better that "conservation". Reducing the runs the cost decreased a lot and the result didn't changed significantly. With respect to crossover and mutation a bigger crossover and lower mutations improved the results. The penalty, as expected, damaged the result, but that doesn't mean the solution will be worse, just that more features choosen are penalised, and it's possible that this will be better than the others regarding to the AUC value.

**Validation of the goodness of the results through the AUC value** 

Function for the AUC with 10 K-fold Cross-Validation to evaluate the attributes given by the filter and wrapper methods, using the models `glm`, `dt`, `nnet`.


```r
crossValidation<-function(all.data, k, formula) {
  
  auc.glm<-list()
  auc.dt<-list()
  auc.nnet<-list()
  cf<-createFolds(all.data$myclass, k=k)
  
  for(i in 1:k) {
    test <-all.data[unlist(cf[i]),]
    train <-all.data[-unlist(cf[i]),]
    train$myclass <-as.factor(train$myclass)
    test$myclass <-as.factor(test$myclass);
    auc.dt<-c(auc.dt,compute.models(train, test, "dt", formula)$auc);
    auc.glm<-c(auc.glm,compute.models(train, test, "glm", formula)$auc);
    auc.nnet<-c(auc.nnet,compute.models(train, test, "nnet", formula)$auc);
  }
  
  #means
  mean.glm=mean(unlist(auc.glm))
  mean.dt=mean(unlist(auc.dt))
  mean.nnet=mean(unlist(auc.nnet))
  means.cv<-list(mean.glm, mean.dt, mean.nnet)
  #boxplot   
  boxplot(as.numeric(auc.glm), as.numeric(auc.dt),as.numeric(auc.nnet), names=c("GLM",  "DT", "NNET"))
  names=c("GLM", "DT", "NNET")
  points(x=1:3,y=means.cv,pch=19,col="green4")
  text(x=1:3,y=means.cv, format(round(as.numeric(means.cv), 2), nsmall = 3), col="mediumseagreen", pos=1)
  
  return(list(glm.mean=mean.glm, dt.mean=mean.dt, nnet.mean=mean.nnet))
}
```

We're going to see the AUC values for the FINAL solutions given by GA algorithm. For extract the list of features corresponding to the solution we're going to use the function implemented in [2.4](#link24)



```r
ga.results<-c(GA.with.penalty,GA.1000.chi, GA.500.chi, GA.500.correlation,GA.500.sbs)
list.appearancesF<-features.frecuency(ga.results)

with.penalty.solution<-list.appearancesF[[1]]
chi1000.solution<-list.appearancesF[[2]]
chi500.solution<-list.appearancesF[[3]]
chi500.correlation<-list.appearancesF[[4]]
chi500.sbs<-list.appearancesF[[5]]

all.data.from.features<-function(train, test, features.in.solution) {
  train.in.solution<-dplyr::select(train, one_of(features.in.solution))
  train.in.solution$myclass<-myclass.train
  test.in.solution<-dplyr::select(test, one_of(features.in.solution))
  test.in.solution$myclass<-myclass.test
  all.data<-rbind(train.in.solution,test.in.solution)
  return(all.data)
}

list.of.features<-list(with.penalty.solution, chi1000.solution, chi500.solution, chi500.correlation, chi500.sbs)
list.of.data<-lapply(list.of.features, function(x) all.data.from.features(ALL.AML_train, ALL.AML_test, x))
list.of.formulas<-lapply(list.of.features, function(x) as.simple.formula(paste(x, collapse="+"), "myclass"))

library(pROC)
```

```
## Type 'citation("pROC")' for a citation.
```

```
## 
## Attaching package: 'pROC'
```

```
## The following objects are masked from 'package:stats':
## 
##     cov, smooth, var
```

```r
set.seed(7)
#With penalty
crossValidation(list.of.data[[1]], 10, list.of.formulas[[1]])
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-80-1.png)<!-- -->

```
## $glm.mean
## [1] 0.9666667
## 
## $dt.mean
## [1] 0.8291667
## 
## $nnet.mean
## [1] 0.5166667
```

```r
#Chi with 100 features
crossValidation(list.of.data[[2]], 10, list.of.formulas[[2]])
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-80-2.png)<!-- -->

```
## $glm.mean
## [1] 0.9416667
## 
## $dt.mean
## [1] 0.9
## 
## $nnet.mean
## [1] 0.5583333
```

```r
#Chi with 500 features
crossValidation(list.of.data[[3]], 10, list.of.formulas[[3]])
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-80-3.png)<!-- -->

```
## $glm.mean
## [1] 0.95
## 
## $dt.mean
## [1] 0.9
## 
## $nnet.mean
## [1] 0.6125
```

```r
#Correlation with 500 features
crossValidation(list.of.data[[4]], 10, list.of.formulas[[4]])
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-80-4.png)<!-- -->

```
## $glm.mean
## [1] 0.9416667
## 
## $dt.mean
## [1] 0.9333333
## 
## $nnet.mean
## [1] 0.5791667
```

```r
#SBS with 500 features
crossValidation(list.of.data[[5]], 10, list.of.formulas[[5]]) 
```

![](P2_FeatureSelection_files/figure-html/unnamed-chunk-80-5.png)<!-- -->

```
## $glm.mean
## [1] 1
## 
## $dt.mean
## [1] 0.8083333
## 
## $nnet.mean
## [1] 0.6083333
```

| Feature selection method | GLM | DT | NNET | 
|--------------------------|----|-----|------|
|Wrapper SBS 500 |  1 | 0.81 | 0.61 |
|Using penalty | 0.97 | 0.83 | 0.52 |
|Filter Chi 500 |  0.95 | 0.9 | 0.61 |
|Filter Correlation 500 | 0.94 | 0.93 | 0.58 |
|Filter Chi 1000 |  0.94 | 0.9 | 0.56 | 



Here we have the real measures of how good and usefuls the selections have been. For *Neural Networks* none of the selections methods give good results cause this particulary one needs large datasets to train the networks, and we have provide a small one. The reuslts for *Linear Regresion* are a sightly better than for *Decision Trees*.

Evidence demostrated: **the top best one is: the Wrapper method SBS**, not only in time and computational cost but also in accuracy.

We have kepp the subset of 100 features with the purpose of proving that the much less features (the half) can give a better result, having also the half of computational cost. Not only it's worse that the one made by the same filter, besides it's the wrost of all: another proof of that **Feature selection is needed**.

And lastly, the result given by the algorithm using **penalty**. Only seeing the fitness value, of course it was the wrost, but 	in effect, as we predicted, it gives more quality features, being the **second best one**.







 
-------




