---
title: 'Final project: Paper replication'
author: "Naira María Chiclana García 44717497T"
date: "February 2019"
output: 
  html_document:
    self_contained: false
    keep_md: true
    toc: true
    toc_float: true
---
#1. About the paper {#paper}

**"OncoLnc: linking TCGA survival data to mRNAs, miRNAs, and lncRNAs"**


**Abstract**: 

*OncoLnc is a tool for interactively exploring survival correlations, and for downloading clinical data coupled to expression data for mRNAs, miRNAs, or long noncoding RNAs (lncRNAs). OncoLnc contains survival data for 8,647 patients from 21 cancer studies performed by The Cancer Genome Atlas (TCGA), along with RNA-SEQ expression for mRNAs and miRNAs from TCGA, and lncRNA expression from MiTranscriptome beta. Storing this data gives users the ability to separate patients by gene expression, and then create publication-quality Kaplan-Meier plots or download the data for further analyses. OncoLnc also stores precomputed survival analyses, allowing users to quickly explore survival correlations for up to 21 cancers in a single click. This resource allows researchers studying a specific gene to quickly investigate if it may have a role in cancer, and the supporting data allows researchers studying a specific cancer to identify the mRNAs, miRNAs, and lncRNAs most correlated with survival, and researchers looking for a novel lncRNA involved with cancer lists of potential candidates. OncoLnc is available at http://www.oncolnc.org. *


Basically, it's a web application where you submit a gene id, and basing on the expression you obtain the value of some parameters (*Cox Coefficient, p-value, FDR, Rank, Median Expresion* and *Mean Expresion*) in addition to a Kaplain Meier survival curve plot.

The used data has been extracted from *TCGA* (The Cancer Genome Atlas). 

The advantages that this program offers are:

- Not having to deal with the complexity and difficulties to acces and fully utilize the data from *TCGA*.

- Focusing on survival analysis with *TCGA* data (altough multiple resources to help researches download TCGA data no one is specific for survival). 

- Allow users to view results from more than one cancer at a time in order to perform a comprehensive survival stdudy with the gene of interest, without missing any interesting correlation. (For not a super extended project we will only do the example with *BRCA*).


**Method:** 

Our objetive is to find a dataset that we can use, of clinical and expression data [2.Data](#data). Then, after some modifications (imputation, grouping levels with *kmeans*..) and extra "experiments" ([3.Clasification with Hierarchical clustering](#hierarchical), we will find the best combination of our clinical data for training a Cox model (checking how good it's by the similarity of it's survival predictions with the real ones) [4.1.Best covariates for Cox Model](#best_parameters_cox). Once that combination of clinical data is found, we will train the Cox model adding all the expression data, and we will compare the results of some genes with the results of the OncoLnc paper tool (given by their respective Cox model)  [4.2. Paper values using expression data](#comp_cox_paper).


#2. Data {#data}

##2.1.BRCA Data provided in paper

We can access to the paper repository in here: https://github.com/OmnesRes/onco_lnc.

As well as in the paper, in the readme, it's said *"This repository already contains the clinical data necessary, but the expression data will need to be downloaded from https://tcga-data.nci.nih.gov/tcga/" * and *"It is also important to note that the names used in TCGA data are outdated" *.

Let's see what contains and how useful it's the clinical data they provide.


```r
setwd("/Users/nairachiclana/Desktop/PFinal-Paper/intento1-oncolog/onco_lnc-master/tcga_data/BRCA/clinical/")
v0<-read.delim("nationwidechildrens.org_clinical_patient_brca.txt", header=T, sep="\t", dec=".")
v1<-read.delim("nationwidechildrens.org_clinical_follow_up_v1.5_brca.txt", header=T, sep="\t", dec=".")
v2<-read.delim("nationwidechildrens.org_clinical_follow_up_v2.1_brca.txt", header=T, sep="\t", dec=".")
v4<-read.delim("nationwidechildrens.org_clinical_follow_up_v4.0_brca.txt", header=T, sep="\t", dec=".")
```

Since there are 4 different files with data, we will join all of them by the colum `bcr_patient_uuid`, the unique identifier of each patient. For that, we will use the `dplyr` function `full_join`.



```r
library(dplyr)
```


```r
v0v1_full<-dplyr::full_join(v0,v1, by = "bcr_patient_uuid")
v0v1v2_full<-dplyr::full_join(v0v1_full,v2, by = "bcr_patient_uuid")
v0v1v2v4_full<-dplyr::full_join(v0v1v2_full,v4, by = "bcr_patient_uuid")
dim(v0v1v2v4_full)
```

```
## [1] 1205  230
```

```r
length(unique(v0v1v2v4_full$bcr_patient_uuid)) 
```

```
## [1] 1099
```

With the full join we have a dataset of 1205 patients and 230 clinical features. Some of the patients (id's) are duplicated.

The covariates are:


```r
names(v0v1v2v4_full)
```

```
##   [1] "bcr_patient_uuid"                                 
##   [2] "bcr_patient_barcode.x"                            
##   [3] "form_completion_date.x"                           
##   [4] "prospective_collection"                           
##   [5] "retrospective_collection"                         
##   [6] "birth_days_to"                                    
##   [7] "gender"                                           
##   [8] "menopause_status"                                 
##   [9] "race"                                             
##  [10] "ethnicity"                                        
##  [11] "history_other_malignancy"                         
##  [12] "history_neoadjuvant_treatment"                    
##  [13] "tumor_status.x"                                   
##  [14] "vital_status.x"                                   
##  [15] "last_contact_days_to.x"                           
##  [16] "death_days_to.x"                                  
##  [17] "radiation_treatment_adjuvant.x"                   
##  [18] "pharmaceutical_tx_adjuvant.x"                     
##  [19] "histologic_diagnosis_other"                       
##  [20] "initial_pathologic_dx_year"                       
##  [21] "age_at_diagnosis"                                 
##  [22] "method_initial_path_dx"                           
##  [23] "method_initial_path_dx_other"                     
##  [24] "surgical_procedure_first"                         
##  [25] "first_surgical_procedure_other"                   
##  [26] "margin_status"                                    
##  [27] "surgery_for_positive_margins"                     
##  [28] "surgery_for_positive_margins_other"               
##  [29] "margin_status_reexcision"                         
##  [30] "axillary_staging_method"                          
##  [31] "axillary_staging_method_other"                    
##  [32] "micromet_detection_by_ihc"                        
##  [33] "lymph_nodes_examined"                             
##  [34] "lymph_nodes_examined_count"                       
##  [35] "lymph_nodes_examined_he_count"                    
##  [36] "lymph_nodes_examined_ihc_count"                   
##  [37] "ajcc_staging_edition"                             
##  [38] "ajcc_tumor_pathologic_pt"                         
##  [39] "ajcc_nodes_pathologic_pn"                         
##  [40] "ajcc_metastasis_pathologic_pm"                    
##  [41] "ajcc_pathologic_tumor_stage"                      
##  [42] "metastasis_site"                                  
##  [43] "metastasis_site_other"                            
##  [44] "er_status_by_ihc.x"                               
##  [45] "er_status_ihc_Percent_Positive.x"                 
##  [46] "er_positivity_scale_used"                         
##  [47] "er_ihc_score.x"                                   
##  [48] "er_positivity_scale_other.x"                      
##  [49] "er_positivity_method.x"                           
##  [50] "pr_status_by_ihc.x"                               
##  [51] "pr_status_ihc_percent_positive.x"                 
##  [52] "pr_positivity_scale_used"                         
##  [53] "pr_positivity_ihc_intensity_score.x"              
##  [54] "pr_positivity_scale_other.x"                      
##  [55] "pr_positivity_define_method.x"                    
##  [56] "her2_status_by_ihc.x"                             
##  [57] "her2_ihc_percent_positive.x"                      
##  [58] "her2_ihc_score.x"                                 
##  [59] "her2_positivity_scale_other.x"                    
##  [60] "her2_positivity_method_text.x"                    
##  [61] "her2_fish_status.x"                               
##  [62] "her2_copy_number.x"                               
##  [63] "cent17_copy_number.x"                             
##  [64] "her2_and_cent17_cells_count.x"                    
##  [65] "her2_cent17_ratio.x"                              
##  [66] "her2_and_cent17_scale_other.x"                    
##  [67] "her2_fish_method.x"                               
##  [68] "new_tumor_event_dx_indicator.x"                   
##  [69] "nte_er_status.x"                                  
##  [70] "nte_er_status_ihc__positive.x"                    
##  [71] "nte_er_ihc_intensity_score.x"                     
##  [72] "nte_er_positivity_other_scale.x"                  
##  [73] "nte_er_positivity_define_method.x"                
##  [74] "nte_pr_status_by_ihc.x"                           
##  [75] "nte_pr_status_ihc__positive.x"                    
##  [76] "nte_pr_ihc_intensity_score.x"                     
##  [77] "nte_pr_positivity_other_scale.x"                  
##  [78] "nte_pr_positivity_define_method.x"                
##  [79] "nte_her2_status.x"                                
##  [80] "nte_her2_status_ihc__positive.x"                  
##  [81] "nte_her2_positivity_ihc_score.x"                  
##  [82] "nte_her2_positivity_other_scale.x"                
##  [83] "nte_her2_positivity_method.x"                     
##  [84] "nte_her2_fish_status.x"                           
##  [85] "nte_her2_signal_number.x"                         
##  [86] "nte_cent_17_signal_number.x"                      
##  [87] "her2_cent17_counted_cells_count.x"                
##  [88] "nte_cent_17_her2_ratio.x"                         
##  [89] "nte_cent17_her2_other_scale.x"                    
##  [90] "nte_her2_fish_define_method.x"                    
##  [91] "anatomic_neoplasm_subdivision"                    
##  [92] "clinical_M"                                       
##  [93] "clinical_N"                                       
##  [94] "clinical_T"                                       
##  [95] "clinical_stage"                                   
##  [96] "days_to_initial_pathologic_diagnosis"             
##  [97] "days_to_patient_progression_free"                 
##  [98] "days_to_tumor_progression"                        
##  [99] "disease_code"                                     
## [100] "extranodal_involvement"                           
## [101] "histological_type"                                
## [102] "icd_10"                                           
## [103] "icd_o_3_histology"                                
## [104] "icd_o_3_site"                                     
## [105] "informed_consent_verified"                        
## [106] "metastatic_tumor_indicator"                       
## [107] "patient_id"                                       
## [108] "project_code"                                     
## [109] "site_of_primary_tumor_other"                      
## [110] "stage_other"                                      
## [111] "tissue_source_site"                               
## [112] "tumor_tissue_site"                                
## [113] "bcr_patient_barcode.y"                            
## [114] "bcr_followup_barcode.x"                           
## [115] "bcr_followup_uuid.x"                              
## [116] "form_completion_date.y"                           
## [117] "radiation_treatment_adjuvant.y"                   
## [118] "tumor_status.y"                                   
## [119] "vital_status.y"                                   
## [120] "last_contact_days_to.y"                           
## [121] "death_days_to.y"                                  
## [122] "new_tumor_event_dx_days_to.x"                     
## [123] "new_tumor_event_radiation_tx.x"                   
## [124] "new_tumor_event_pharmaceutical_tx.x"              
## [125] "nte_er_status.y"                                  
## [126] "nte_er_status_ihc__positive.y"                    
## [127] "nte_er_ihc_intensity_score.y"                     
## [128] "nte_er_positivity_other_scale.y"                  
## [129] "nte_er_positivity_define_method.y"                
## [130] "nte_pr_status_by_ihc.y"                           
## [131] "nte_pr_status_ihc__positive.y"                    
## [132] "nte_pr_ihc_intensity_score.y"                     
## [133] "nte_pr_positivity_other_scale.y"                  
## [134] "nte_pr_positivity_define_method.y"                
## [135] "nte_her2_status.y"                                
## [136] "nte_her2_status_ihc__positive.y"                  
## [137] "nte_her2_positivity_ihc_score.y"                  
## [138] "nte_her2_positivity_other_scale.y"                
## [139] "nte_her2_positivity_method.y"                     
## [140] "nte_her2_fish_status.y"                           
## [141] "nte_her2_signal_number.y"                         
## [142] "nte_cent_17_signal_number.y"                      
## [143] "her2_cent17_counted_cells_count.y"                
## [144] "nte_cent_17_her2_ratio.y"                         
## [145] "nte_cent17_her2_other_scale.y"                    
## [146] "nte_her2_fish_define_method.y"                    
## [147] "cent17_copy_number.y"                             
## [148] "days_to_additional_surgery_locoregional_procedure"
## [149] "days_to_additional_surgery_metastatic_procedure.x"
## [150] "days_to_last_known_alive"                         
## [151] "er_ihc_score.y"                                   
## [152] "er_positivity_method.y"                           
## [153] "er_positivity_scale_other.y"                      
## [154] "er_status_by_ihc.y"                               
## [155] "er_status_ihc_Percent_Positive.y"                 
## [156] "her2_and_cent17_cells_count.y"                    
## [157] "her2_and_cent17_scale_other.y"                    
## [158] "her2_cent17_ratio.y"                              
## [159] "her2_copy_number.y"                               
## [160] "her2_fish_method.y"                               
## [161] "her2_fish_status.y"                               
## [162] "her2_ihc_percent_positive.y"                      
## [163] "her2_ihc_score.y"                                 
## [164] "her2_positivity_method_text.y"                    
## [165] "her2_positivity_scale_other.y"                    
## [166] "her2_status_by_ihc.y"                             
## [167] "new_tumor_event_surgery.x"                        
## [168] "new_tumor_event_surgery_met"                      
## [169] "pr_positivity_define_method.y"                    
## [170] "pr_positivity_ihc_intensity_score.y"              
## [171] "pr_positivity_scale_other.y"                      
## [172] "pr_status_by_ihc.y"                               
## [173] "pr_status_ihc_percent_positive.y"                 
## [174] "targeted_molecular_therapy.x"                     
## [175] "bcr_patient_barcode.x.x"                          
## [176] "bcr_followup_barcode.y"                           
## [177] "bcr_followup_uuid.y"                              
## [178] "form_completion_date.x.x"                         
## [179] "radiation_treatment_adjuvant.x.x"                 
## [180] "tumor_status.x.x"                                 
## [181] "vital_status.x.x"                                 
## [182] "last_contact_days_to.x.x"                         
## [183] "death_days_to.x.x"                                
## [184] "new_tumor_event_dx_indicator.y"                   
## [185] "new_tumor_event_type"                             
## [186] "new_tumor_event_site"                             
## [187] "new_tumor_event_site_other"                       
## [188] "new_tumor_event_dx_days_to.y"                     
## [189] "new_tumor_event_radiation_tx.y"                   
## [190] "new_tumor_event_pharmaceutical_tx.y"              
## [191] "nte_er_status"                                    
## [192] "nte_er_status_ihc__positive"                      
## [193] "nte_er_positivity_scale_used"                     
## [194] "nte_er_ihc_intensity_score"                       
## [195] "nte_er_positivity_other_scale"                    
## [196] "nte_er_positivity_define_method"                  
## [197] "nte_pr_status_by_ihc"                             
## [198] "nte_pr_status_ihc__positive"                      
## [199] "nte_pr_positivity_scale_used"                     
## [200] "nte_pr_ihc_intensity_score"                       
## [201] "nte_pr_positivity_other_scale"                    
## [202] "nte_pr_positivity_define_method"                  
## [203] "nte_her2_status"                                  
## [204] "nte_her2_status_ihc__positive"                    
## [205] "nte_her2_positivity_ihc_score"                    
## [206] "nte_her2_positivity_other_scale"                  
## [207] "nte_her2_positivity_method"                       
## [208] "nte_her2_fish_status"                             
## [209] "nte_her2_signal_number"                           
## [210] "nte_cent_17_signal_number"                        
## [211] "her2_cent17_counted_cells_count"                  
## [212] "nte_cent_17_her2_ratio"                           
## [213] "nte_cent17_her2_other_scale"                      
## [214] "nte_her2_fish_define_method"                      
## [215] "days_to_additional_surgery_metastatic_procedure.y"
## [216] "followup_reason"                                  
## [217] "new_tumor_event_surgery.y"                        
## [218] "targeted_molecular_therapy.y"                     
## [219] "bcr_patient_barcode.y.y"                          
## [220] "bcr_followup_barcode"                             
## [221] "bcr_followup_uuid"                                
## [222] "form_completion_date.y.y"                         
## [223] "followup_lost_to"                                 
## [224] "radiation_treatment_adjuvant.y.y"                 
## [225] "pharmaceutical_tx_adjuvant.y"                     
## [226] "tumor_status.y.y"                                 
## [227] "vital_status.y.y"                                 
## [228] "last_contact_days_to.y.y"                         
## [229] "death_days_to.y.y"                                
## [230] "new_tumor_event_dx_indicator"
```

```r
full.df<-v0v1v2v4_full
```



###2.1.1.Cox using paper raw data {#paper_data}


First of all, before any modifications, we will try to use the raw data for Cox.

It's	unavoidable at least,  correct the covariates *living status* and *time*.

We will remove the level  *"CDE_ID:5"* which does not mean anything to us, and transform the status level to 0 (censured) and 1 (death).


```r
#Covariate living status
full.df<-full.df[2:nrow(full.df),] #first row are repeated names of covariates
levels(full.df$vital_status.x)
```

```
## [1] "Alive"        "CDE_ID:5"     "Dead"         "vital_status"
```

```r
full.df<-full.df[full.df$vital_status.x!="CDE_ID:5",]
full.df$vital_status.x<-as.numeric(factor(full.df$vital_status.x, labels=c("0", "1")))  -1
str(full.df$vital_status.x)
```

```
##  num [1:1203] 0 0 0 0 0 0 0 0 0 0 ...
```

```r
table(full.df$vital_status.x)
```

```
## 
##    0    1 
## 1097  106
```

We have 1097 censured patients and 106 patients whose event (death) have occured.

After removing the not avaliable rows for time, we will also parse the time in days to numeric so we can use it in survival functions.


```r
#Covariate time
full.df<-full.df[full.df$last_contact_days_to.x!="[Not Available]",]
full.df$last_contact_days_to.x<-as.numeric(full.df$last_contact_days_to.x)
str(full.df$last_contact_days_to.x)
```

```
##  num [1:1097] 408 408 412 66 66 62 62 156 133 132 ...
```

**Cox model** using all parameters:


```r
library("survival")
library("survminer")
```

```r
full.cox<-coxph(Surv(last_contact_days_to.x, vital_status.x)~., data=full.df)
full.cox$loglik
```

```
## [1] 0 0
```

All the values are *<NA>* and the likelihood is 0. We can't do Cox Model with this data. This is probably because the majority of the data is empty or *Not avaliable*.

Before discarding it, we will see how the survival curve looks.


```r
ggsurvplot(survfit(Surv(last_contact_days_to.x,vital_status.x)~1, data=full.df, type="kaplan-meier"))
```

![](paper_files/figure-html/unnamed-chunk-9-1.png)<!-- -->


Definitely, we can't do anything with the data in this way. 

We will try again but **removing the rows with missing values**:

First, we will make a list with the names of covariates whose have any missing value, and then we will iterate over these columns removing the rows with any missing data.


```r
#Covariates which contain not avaliable values
variables_need_imputation<-list()
for(variable in colnames(full.df)) {
  if("[Not Available]" %in% levels(full.df[,variable])) {
    variables_need_imputation<-c(variables_need_imputation, variable)
  }
}
variables_need_imputation<-unlist(variables_need_imputation) 

#Iterate over the columns which need imputation and remove the rows with missing values
full.df.imp<-full.df
for(variable in variables_need_imputation) {
  col<-full.df.imp[,variable]
  full.df.imp<-full.df.imp[col!="[Not Available]",]
}

dim(full.df)
```

```
## [1] 1097  230
```

```r
dim(full.df.imp)
```

```
## [1]   0 230
```

While the original dataset had 1097 rows, the imputed dataset has no one row remaining. This is completely **unusable**.

Like it was said in the paper, we have checked that this data is outdated and useless. We will get the data directly from *TCGA* ourserlves.


##2.2. BRCA Data from TCGA


Using the package `TCGAretriever` we will search the clinical data we want from BRCA patients.



```r
library(TCGAretriever) 
```



```r
all_studies<-get_cancer_studies() #identifier, tittle and description of the study
head(all_studies[,2],20) #tittles of first 20 studies
```

```
##  [1] "Acinar Cell Carcinoma of the Pancreas (Johns Hopkins, J Pathol 2014)"
##  [2] "Acute Lymphoblastic Leukemia (St Jude, Nat Genet 2015)"              
##  [3] "Acute Lymphoblastic Leukemia (St Jude, Nat Genet 2016)"              
##  [4] "Acute Myeloid Leukemia (TCGA, NEJM 2013)"                            
##  [5] "Acute Myeloid Leukemia (TCGA, PanCancer Atlas)"                      
##  [6] "Acute Myeloid Leukemia (TCGA, Provisional)"                          
##  [7] "Adenoid Cystic Carcinoma (FMI, Am J Surg Pathl. 2014)"               
##  [8] "Adenoid Cystic Carcinoma (MDA, Clin Cancer Res 2015)"                
##  [9] "Adenoid Cystic Carcinoma (MGH, Nat Gen 2016)"                        
## [10] "Adenoid Cystic Carcinoma (MSKCC, Nat Genet 2013)"                    
## [11] "Adenoid Cystic Carcinoma (Sanger/MDA, JCI 2013)"                     
## [12] "Adenoid Cystic Carcinoma of the Breast (MSKCC, J Pathol. 2015)"      
## [13] "Adrenocortical Carcinoma (TCGA, PanCancer Atlas)"                    
## [14] "Adrenocortical Carcinoma (TCGA, Provisional)"                        
## [15] "Ampullary Carcinoma (Baylor College of Medicine, Cell Reports 2016)" 
## [16] "Bladder Cancer (MSKCC, Eur Urol 2014)"                               
## [17] "Bladder Cancer (MSKCC, JCO 2013)"                                    
## [18] "Bladder Cancer (MSKCC, Nat Genet 2016)"                              
## [19] "Bladder Cancer (TCGA, Cell 2017)"                                    
## [20] "Bladder Urothelial Carcinoma (BGI, Nat Genet 2013)"
```

```r
#filter rows from all_studies which name contains breast
df.breast<-dplyr::select(filter(all_studies, grepl('Breast', all_studies$name)),everything()) 
df.breast$cancer_study_id
```

```
##  [1] "acbc_mskcc_2015"              "brca_metabric"               
##  [3] "breast_msk_2018"              "bfn_duke_nus_2015"           
##  [5] "brca_bccrc"                   "brca_broad"                  
##  [7] "brca_sanger"                  "brca_tcga_pub2015"           
##  [9] "brca_tcga_pub"                "brca_tcga_pan_can_atlas_2018"
## [11] "brca_tcga"                    "brca_bccrc_xenograft_2014"   
## [13] "brca_mbcproject_wagle_2017"
```

```r
#from all the studies, we choose brca_tcga->brca_tcga_all
all_case_lists<-get_case_lists("brca_tcga")
colnames(all_case_lists)
```

```
## [1] "case_list_id"          "case_list_name"        "case_list_description"
## [4] "cancer_study_id"       "case_ids"
```

```r
all_case_lists[,1:2] #id and name
```

```
##                        case_list_id
## 1           brca_tcga_3way_complete
## 2               brca_tcga_sequenced
## 3         brca_tcga_methylation_all
## 4                     brca_tcga_all
## 5  brca_tcga_protein_quantification
## 6                     brca_tcga_cna
## 7        brca_tcga_methylation_hm27
## 8       brca_tcga_methylation_hm450
## 9                    brca_tcga_mrna
## 10        brca_tcga_rna_seq_v2_mrna
## 11                   brca_tcga_rppa
## 12                 brca_tcga_cnaseq
##                                       case_list_name
## 1                                All Complete Tumors
## 2                               All Sequenced Tumors
## 3            All tumor samples with methylation data
## 4                                         All Tumors
## 5                 Protein Quantification (Mass Spec)
## 6                        Tumor Samples with CNA data
## 7         Tumor Samples with methylation data (HM27)
## 8        Tumor Samples with methylation data (HM450)
## 9  Tumor Samples with mRNA data (Agilent microarray)
## 10         Tumor Samples with mRNA data (RNA Seq V2)
## 11                      Tumor Samples with RPPA data
## 12        Tumor Samples with sequencing and CNA data
```

```r
df.tcga<-get_clinical_data("brca_tcga_all")
```

```
## Warning in (function (..., deparse.level = 1) : number of columns of result
## is not a multiple of vector length (arg 674)
```

After a preview of all the avaliable covariates, we will choose those which seem more relevants for our study.


```r
str(df.tcga$AGE) #age
```

```
##  chr [1:1105] "62" "59" "70" "42" "69" "61" "50" "80" "59" "65" "79" ...
```

```r
str(df.tcga$MENOPAUSE_STATUS) #menostat
```

```
##  chr [1:1105] "Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy)" ...
```

```r
str(df.tcga$AJCC_PATHOLOGIC_TUMOR_STAGE) #tgrade
```

```
##  chr [1:1105] "Stage IIA" "Stage IIIC" "Stage IIIA" "Stage IIB" ...
```

```r
str(df.tcga$LYMPH_NODE_EXAMINED_COUNT) #pnodes
```

```
##  chr [1:1105] "2" "15" "17" "8" "4" "4" "30" "10" "" "4" "2" "10" "12" ...
```

```r
str(df.tcga$ER_STATUS_BY_IHC)  #strec
```

```
##  chr [1:1105] "Positive" "Positive" "Positive" "Positive" "Positive" ...
```

```r
str(df.tcga$PR_STATUS_BY_IHC) #progrec
```

```
##  chr [1:1105] "Positive" "Negative" "Positive" "Positive" "Positive" ...
```

```r
str(df.tcga$IHC_HER2) #her2
```

```
##  chr [1:1105] "Negative" "" "Positive" "Positive" "Equivocal" ...
```

```r
str(df.tcga$OS_MONTHS) #time (months)
```

```
##  chr [1:1105] "10.28" "26.02" "5.65" "40.47" "24.47" "97.77" "85.58" ...
```

```r
str(df.tcga$OS_STATUS) #cens
```

```
##  chr [1:1105] "LIVING" "LIVING" "DECEASED" "LIVING" "LIVING" "LIVING" ...
```

```r
str(df.tcga$CASE_ID) #patient id
```

```
##  chr [1:1105] "TCGA-A7-A3J0-01" "TCGA-OL-A66N-01" "TCGA-AQ-A0Y5-01" ...
```



###2.2.1. Preparation of the data

We have seen in [Cox using paper raw data](#paper_data) how much having missing values can affect, so we will remove them (the row to where they belong) from our new dataset.


#### Removing missing values: 

We will remove the rows with empty or indetermiate values.


```r
#Remove inteterminate and empty values 
  #remove os_months=0
df.tcga<-df.tcga[df.tcga$OS_MONTHS>0,]
  #remove menostat with values Indeterminate or empty
df.tcga<-df.tcga[as.factor(df.tcga$MENOPAUSE_STATUS)!="Indeterminate (neither Pre or Postmenopausal)",]
df.tcga<-df.tcga[as.factor(df.tcga$MENOPAUSE_STATUS)!="",]
  #remove tgrade empty
df.tcga<-df.tcga[as.factor(df.tcga$AJCC_PATHOLOGIC_TUMOR_STAGE)!="",]
  #remove HER, PR, ER empty or indeterminate
df.tcga<-df.tcga[as.factor(df.tcga$ER_STATUS_BY_IHC)!="",]
df.tcga<-df.tcga[as.factor(df.tcga$PR_STATUS_BY_IHC)!="",]
df.tcga<-df.tcga[as.factor(df.tcga$PR_STATUS_BY_IHC)!="Indeterminate",]
df.tcga<-df.tcga[as.factor(df.tcga$IHC_HER2)!="",]
df.tcga<-df.tcga[as.factor(df.tcga$IHC_HER2)!="Indeterminate",]

#Create dataframe
df.tcga$OS_MONTHS<-as.numeric(df.tcga$OS_MONTHS)
df.final<-data.frame(df.tcga$CASE_ID, df.tcga$AGE,df.tcga$MENOPAUSE_STATUS,df.tcga$AJCC_PATHOLOGIC_TUMOR_STAGE,df.tcga$LYMPH_NODES_EXAMINED_HE_COUNT,df.tcga$ER_STATUS_BY_IHC,df.tcga$PR_STATUS_BY_IHC,df.tcga$IHC_HER2,df.tcga$OS_MONTHS,df.tcga$OS_STATUS)
#name dataframe columns
names(df.final)<-c("patient_id","age", "menostat",  "tgrade", "pnodes", "ER", "PR", "HER2", "time", "cens")
dim(df.final)
```

```
## [1] 817  10
```

After removing the missing and indeterminate values, we still have 817 patients remaining, an acceptable quantity to work with.

#### Imputation of numerical and categorical covariates: 

For numerical, we will replace the missing values of number of nodes and age with their mean.
For categoricals, we will replace them with the most suitable value depending of each one.


```r
impute.cols<-function(df) {
  
  #NUMERIC->MEAN
  imputar.numericos<-function(df, col) {
    mean.value<-mean(as.integer(na.omit(df[,col])))
    df[,col][df[,col]==""]<-floor(mean.value)
    df[,col]<-as.numeric(as.character(df[,col])) 
    return(df[,col])
  }

  df$age<-imputar.numericos(df,"age")
  df$pnodes<-imputar.numericos(df,"pnodes")
  
  #CATEROTICALS
  cat("LEVELS OF CATEGORICALS: \n")
  #menostat
  mode.menostat<-names(table(df$menostat)[which.max(table(df$menostat))])
  mode.menostat
  cat("Menostat levels: ", levels(df$menostat), "\n") # Peri as (Post)
  df$menostat<-factor(df$menostat, labels=c("Post", "Post", "Pre"))
  #tgrade
  mode.tgrade<-names(table(df$tgrade)[which.max(table(df$tgrade))])
  mode.tgrade
  cat("Tumor grade levels: ", levels(df$tgrade), "\n") #X as mode (II)
  df$tgrade<-factor(df$tgrade, labels=c("I","I", "I", "II", "II", "II", "III", "III", "III", "III", "IV", "II"))
  #cens 
  cat("cens levels: ", levels(df.final$cens), "\n")
  df$cens<-as.numeric(factor(df$cens, labels=c("1", "0"))) -1 #factor starts at 0
  #ER
  cat("Estrogen Receptor levels: ", levels(df$ER), "\n") 
  df$ER<-factor(df$ER, labels=c("Negative", "Positive"))
  #HER2
  cat("HER2 levels: ",levels(df$HER2), "\n") #Equivocal as Negative
  df$HER2<-factor(df$HER2, labels=c("Negative", "Negative", "Positive"))
  #PR
  cat("Progesterone levels: ", levels(df$PR), "\n")

df$PR<-factor(df$PR, labels=c("Negative", "Positive"))
  
  return(df)
}

df.final.imp<-impute.cols(df.final)
```

```
## LEVELS OF CATEGORICALS: 
## Menostat levels:  Peri (6-12 months since last menstrual period) Post (prior bilateral ovariectomy OR >12 mo since LMP with no prior hysterectomy) Pre (<6 months since LMP AND no prior bilateral ovariectomy AND not on estrogen replacement) 
## Tumor grade levels:  Stage I Stage IA Stage IB Stage II Stage IIA Stage IIB Stage III Stage IIIA Stage IIIB Stage IIIC Stage IV Stage X 
## cens levels:  DECEASED LIVING 
## Estrogen Receptor levels:  Negative Positive 
## HER2 levels:  Equivocal Negative Positive 
## Progesterone levels:  Negative Positive
```

```r
df.final.imp[1:20,]
```

```
##         patient_id age menostat tgrade pnodes       ER       PR     HER2
## 1  TCGA-A7-A3J0-01  62     Post     II      0 Positive Positive Negative
## 2  TCGA-AQ-A0Y5-01  70     Post    III      5 Positive Positive Positive
## 3  TCGA-E9-A22H-01  42      Pre     II      1 Positive Positive Positive
## 4  TCGA-BH-A0EB-01  69     Post      I      0 Positive Positive Negative
## 5  TCGA-AR-A24S-01  61     Post      I      0 Positive Positive Negative
## 6  TCGA-AR-A254-01  50      Pre    III      9 Positive Positive Positive
## 7  TCGA-EW-A1IW-01  80     Post     II      2 Positive Positive Positive
## 8  TCGA-A8-A08F-01  59     Post    III      6 Positive Positive Negative
## 9  TCGA-AR-A24H-01  65     Post     II      0 Positive Positive Negative
## 10 TCGA-LL-A5YP-01  49     Post     II      2 Positive Negative Negative
## 11 TCGA-A7-A3IZ-01  62     Post     II      6 Positive Negative Negative
## 12 TCGA-AO-A0JC-01  64     Post     II      0 Positive Positive Negative
## 13 TCGA-E2-A15R-01  64     Post     II      3 Positive Positive Negative
## 14 TCGA-E9-A1N9-01  58     Post     II      0 Negative Positive Positive
## 15 TCGA-A8-A07Z-01  85     Post     II      6 Positive Positive Negative
## 16 TCGA-BH-A0E7-01  79     Post     II      1 Positive Positive Negative
## 17 TCGA-BH-A0DV-01  54     Post    III      9 Positive Positive Negative
## 18 TCGA-D8-A27N-01  36      Pre    III      4 Positive Positive Positive
## 19 TCGA-BH-A0B3-01  53     Post     II      2 Negative Negative Negative
## 20 TCGA-BH-A0DK-01  49      Pre     II      0 Positive Positive Negative
##      time cens
## 1   10.28    1
## 2    5.65    0
## 3   40.47    1
## 4   24.47    1
## 5   97.77    1
## 6   85.58    1
## 7   12.19    1
## 8   32.98    1
## 9  160.78    1
## 10  14.78    1
## 11  10.58    1
## 12  50.82    1
## 13  56.90    1
## 14  36.17    1
## 15  45.04    1
## 16  44.78    1
## 17  67.81    1
## 18  17.05    1
## 19  39.52    1
## 20  13.90    1
```

```r
BC<-df.final.imp
```

#### Stratification:

**Groups defined in medical literature:**

Luminal:

- LuminalA: ER positive and HER2 negative (30-70% of patients).

- Luminal B: ER positive. (10-20% of patients).

- TN (Triple negative/basal-like): ER-negative and PR-negative. (15-20% of patients).


```r
#LUMINAL
make.luminal<-function(df) {
  vectorluminals <-list()
  for(i in 1:nrow(df)) {
    hormone_receptor_positive=(df[i,]$ER=="Positive" || df[i,]$PR=="Positive")
    hormone_receptor_negative=(df[i,]$ER=="Negative" && df[i,]$PR=="Negative")
    if(hormone_receptor_positive && df[i,]$HER2=="Negative") vectorluminals[i]<-"LuminalA"
    else if(hormone_receptor_positive)  vectorluminals[i]<-"LuminalB" 
    else if(df[i,]$ER=="Negative" && df[i,]$PR=="Negative" && df[i,]$HER2=="Negative") vectorluminals[i]<-"TN"
    else if(hormone_receptor_negative && df[i,]$HER2=="Positive") vectorluminals[i]<-"HER2-enriched"
    else if(hormone_receptor_positive && df[i,]$HER2=="Negative") vectorluminals[i]<-"Normal-line"
  }
  df$Luminal<-as.factor(unlist(vectorluminals))
  return(df)
}

BC<-make.luminal(BC)
prop.table(table(BC$Luminal))
```

```
## 
## HER2-enriched      LuminalA      LuminalB            TN 
##    0.04161567    0.65850673    0.13463892    0.16523868
```

The stratification of Luminal agree with the description and propotions expected.


```r
#NODAL AND NODE STATUS
make.nodalStatus<-function(df) {
#five groups: 0, negative; 1, one positive; 2, two to three positive; 3, four to nine positive; and 4, >9 nine positive
  vectorGroups<-list()
  nodeSt<-list()
  for(n in 1:length(df$pnodes)) {
    numnod<-df$pnodes[n]
    if(numnod==0) {
      vectorGroups[n]="G0"
      nodeSt[n]="Node negative"
    } 
    else if(numnod>0) {
      nodeSt[n]="Node positive"
      if(numnod==1) vectorGroups[n]="G1"
      else if(numnod==2 || numnod==3) vectorGroups[n]="G2"
      else if(numnod>=4 &&numnod<=9) vectorGroups[n]="G3"
      else if(numnod>9) vectorGroups[n]="G4"
    }
  }
  df$NodalStatus<-as.factor(unlist(vectorGroups))
  df$NodeStatus<-as.factor(unlist(nodeSt))
  return(df)
}

BC<-make.nodalStatus(BC)

#Otra agrupación distinta
max=max(BC$pnodes)
BC$NodalStatus2<-cut(BC$pnodes, breaks=c(0,1,3,5, max), include.lowest=T)
```


**Not known stratifications with kmeans: age**

Like we saw in the *Lab Project 1: Clustering methods to sporulation yeast problem* we can choose between three different clustering tipes (average, complete and single) and compute the *Silhouette index* as a function of the number of groups. The *Silhouette index* will give us an idea of the coesion and separation of our clusters, if it's well matched.

We will try to find the best number of clusters of the age covariate in our dataset.


```r
library(devtools)
library(cluster)
library(NbClust)
library(factoextra)
library(ggplot2)
library(gridExtra)
library(FunCluster)
```


```r
optimal.number.of.clusters<-function(data, covariate, clustering.types) {
  
  data$d1<-rep(1,length(data[,covariate]))
  data.num=model.matrix(d1+age~age, data)
  
  d<-dist(data, method="euclidean") 
  for(m in 1:length(clustering.types)) {assign(paste("p",m,sep=""),fviz_nbclust(data.frame(data[,covariate], data$d1), FUNcluster=hcut, method ="silhouette", diss=d, hc_method=clustering.types[m])) }
  grid.arrange(p1, p2, p3, ncol = 2)
}

clustering.types<-c("average", "complete", "single")
optimal.number.of.clusters(BC, "age", clustering.types)
```

![](paper_files/figure-html/unnamed-chunk-19-1.png)<!-- -->


The optimal k is 2 for all the three measures. The worst result is for the measure *"single"*.

Knowing this, we will make two clusters for age using **k-means**:


```r
data<-BC
data$d1<-rep(1,length(data$age))
data.num=model.matrix(d1+age~age, data)
data$cl= kmeans(data.num,2)$cluster
BC$ClusteredAge<-cut(BC$age, breaks=c(min(data$age[data$cl==1]),min(BC$age[data$cl==2]),max(data$age[data$cl==2])))
table(BC$ClusteredAge)
```

```
## 
## (26,60] (60,90] 
##     449     367
```

The distribution of data is quite similar, that is good.



```r
head(BC)
```

```
##        patient_id age menostat tgrade pnodes       ER       PR     HER2
## 1 TCGA-A7-A3J0-01  62     Post     II      0 Positive Positive Negative
## 2 TCGA-AQ-A0Y5-01  70     Post    III      5 Positive Positive Positive
## 3 TCGA-E9-A22H-01  42      Pre     II      1 Positive Positive Positive
## 4 TCGA-BH-A0EB-01  69     Post      I      0 Positive Positive Negative
## 5 TCGA-AR-A24S-01  61     Post      I      0 Positive Positive Negative
## 6 TCGA-AR-A254-01  50      Pre    III      9 Positive Positive Positive
##    time cens  Luminal NodalStatus    NodeStatus NodalStatus2 ClusteredAge
## 1 10.28    1 LuminalA          G0 Node negative        [0,1]      (60,90]
## 2  5.65    0 LuminalB          G3 Node positive        (3,5]      (60,90]
## 3 40.47    1 LuminalB          G1 Node positive        [0,1]      (26,60]
## 4 24.47    1 LuminalA          G0 Node negative        [0,1]      (60,90]
## 5 97.77    1 LuminalA          G0 Node negative        [0,1]      (60,90]
## 6 85.58    1 LuminalB          G3 Node positive       (5,29]      (26,60]
```

```r
BC.no.cens<-BC[ , !(names(BC) %in% "cens")]
```

## 2.3. Expression Data {#expression_data}

The paper uses three types of expression data: *mRNAs, miRNAs and long noncoding RNAs (lncRNAs)*. We will use *mRNA* expression data of genes for $BRCA$ obtained from TCGA.


```r
library(readxl)
library(RTCGAToolbox)
library(tibble)
library(dplyr)
library(limma)
library(edgeR)
library(calibrate)

brcaData <- getFirehoseData(dataset="BRCA", runDate="20160128",gistic2Date="20160128",forceDownload=F, clinical =TRUE, RNASeq2GeneNorm  =TRUE)
brca_rnaseq <- getData(brcaData,type = "RNASeq2GeneNorm")
brca_rnaseq.tumour <- brca_rnaseq[, which(as.numeric(substr(colnames(brca_rnaseq), 14,15)) < 10)]
colnames(brca_rnaseq.tumour) <- substr(colnames(brca_rnaseq.tumour), 1,12)
brca_rnaseq.tumour <- brca_rnaseq.tumour[, !duplicated(colnames(brca_rnaseq.tumour))]
temp <- tempfile()
download.file("https://media.nature.com/original/nature-assets/nature/journal/v490/n7418/extref/nature11412-s2.zip",temp)
unzip(temp)
sample_data <- read_excel("nature11412-s2/Supplementary Tables 1-4.xls", sheet = 1, skip = 1)
unlink("nature11412-s2",recursive = T,force = T)
rm(temp)
unlink("20160128*",force = T)

tnbc_samples <- sample_data %>% dplyr::filter(`ER Status` == "Negative" & `PR Status` == "Negative" & `HER2 Final Status` == "Negative" & `PAM50 mRNA` != "Luminal A")
tnbc_barcodes <- tnbc_samples$`Complete TCGA ID`
luminal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Luminal A")
luminal_barcodes <- luminal_samples$`Complete TCGA ID`
basal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Basal-like")
basal_barcodes <- basal_samples$`Complete TCGA ID`

brca_rnaseq.tnbc <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% tnbc_barcodes)]
brca_rnaseq.luminal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% luminal_barcodes)]
brca_rnaseq.basal <- brca_rnaseq.tumour[, which(colnames(brca_rnaseq.tumour) %in% basal_barcodes)]

d1 = brca_rnaseq.tnbc
d2 = brca_rnaseq.luminal

rnaseq.for.de <- cbind(d1, d2)
counts = rnaseq.for.de[apply(rnaseq.for.de,1,function(x) sum(x==0))<ncol(rnaseq.for.de)*0.8,]

df.l <- data_frame("sample" = colnames(d1), "status" = rep(0, length(colnames(d1))) )
df.t <- data_frame("sample" = colnames(d2), "status" = rep(1, length(colnames(d2))) )
df <- rbind(df.t,df.l)
design <- model.matrix(~ status, data = df)

dge <- DGEList(counts=counts)
A <- rowSums(dge$counts)
isexpr <- A > 100 # Keeping genes with total counts more than 100.
dge <- calcNormFactors(dge)
v <- voom(dge[isexpr,], design, plot=FALSE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

diff.exp.df <- topTable(fit, coef = "status", n = Inf, sort = "p", p = 0.01) # Positive log-fold-changes mean higher expression in d1
diff.exp.df$gene.name <- rownames(diff.exp.df)
```


Once we have the expression data, we will take only the overexpressed one. For difference that, we have defined  a $Log For Change$ upper or lower of $2$ and a $p-value$ less than $0.01$.


```r
menor <- diff.exp.df[diff.exp.df$logFC<(-2),]
mayor <- diff.exp.df[diff.exp.df$logFC>2,]
sobreexp <- rbind(menor,mayor)
sobreexp <- sobreexp[ sobreexp$adj.P.Val<0.01,]
genes.model<- data.frame(t(counts[sobreexp$gene.name,]))
```

Having that data, we will make a dataset `BC.expression` which contains the clinical data we already have and the expression data of the genes whose coincidence is found  in the patient id.


```r
renamed.id<-vector()
for(id in 1:length(BC$patient_id)) {
  renamed.id[id]<-paste(strsplit(as.character(BC$patient_id)[id], "-")[[1]][1:3], collapse="-")
}

BC$`Complete TCGA ID`<- renamed.id

genes.model$`Complete TCGA ID` <- rownames(genes.model)
compact<- sample_data[c(sample_data$`Complete TCGA ID`) %in% genes.model$`Complete TCGA ID`,]
compact<- BC[c(BC$`Complete TCGA ID`) %in% genes.model$`Complete TCGA ID`,]
join <- inner_join(compact, genes.model)

clinicos.completos <- sample_data
clinicos.reducidos <- compact
BC.expression <- join
```





```r
dim(BC)
```

```
## [1] 817  15
```

```r
dim(BC.expression)
```

```
## [1] 270 142
```

```r
#Clinical data and genes
names(BC.expression)
```

```
##   [1] "patient_id"       "age"              "menostat"        
##   [4] "tgrade"           "pnodes"           "ER"              
##   [7] "PR"               "HER2"             "time"            
##  [10] "cens"             "Luminal"          "NodalStatus"     
##  [13] "NodeStatus"       "NodalStatus2"     "ClusteredAge"    
##  [16] "Complete TCGA ID" "AGR3"             "DNALI1"          
##  [19] "CAPN8"            "C6orf97"          "GPR77"           
##  [22] "SERPINA11"        "FSIP1"            "TFF3"            
##  [25] "CACNA2D2"         "SYTL5"            "FLJ45983"        
##  [28] "TFF1"             "PRR15"            "LOC145837"       
##  [31] "DACH1"            "C1orf64"          "KCNK15"          
##  [34] "NAT1"             "ANKRD30B"         "MAPT"            
##  [37] "ERBB4"            "ESR1"             "SYT9"            
##  [40] "TSPAN1"           "NEK10"            "ZMYND10"         
##  [43] "DNAJC12"          "ABCC8"            "BCAS1"           
##  [46] "SLC44A4"          "CT62"             "AGR2"            
##  [49] "GP2"              "CHAD"             "PGR"             
##  [52] "GRIK3"            "DEGS2"            "SCUBE2"          
##  [55] "LOC100130148"     "AKR7A3"           "TTC36"           
##  [58] "PTPRT"            "ADAMTS15"         "GATA3"           
##  [61] "IGFALS"           "GREB1"            "ANKRD30A"        
##  [64] "CA12"             "FOXA1"            "MLPH"            
##  [67] "TMPRSS6"          "AFF3"             "F7"              
##  [70] "LOC100128977"     "C20orf114"        "PGLYRP2"         
##  [73] "TPRG1"            "RGS22"            "GRPR"            
##  [76] "AR"               "SPDEF"            "CYP2B7P1"        
##  [79] "SERPINA5"         "LOC728606"        "GFRA1"           
##  [82] "EEF1A2"           "SLC7A2"           "TUBA3E"          
##  [85] "LRG1"             "CYP4B1"           "CYP4Z2P"         
##  [88] "KCNC2"            "MS4A8B"           "SYT1"            
##  [91] "EMX1"             "HMGCS2"           "C1orf173"        
##  [94] "LOC389033"        "SERPINA6"         "TMC5"            
##  [97] "HEPACAM2"         "RIMS4"            "PHGR1"           
## [100] "CHST8"            "CEACAM6"          "CEACAM5"         
## [103] "ABCC11"           "SLC30A8"          "BPIL1"           
## [106] "CPB1"             "KCNJ3"            "C10orf82"        
## [109] "CDC20B"           "SCGB2A2"          "DHRS2"           
## [112] "ABCC13"           "BMPR1B"           "PIP"             
## [115] "SCGB1D2"          "VSTM2A"           "CST9"            
## [118] "CLEC3A"           "FZD9"             "ART3"            
## [121] "CXorf49B"         "CXorf61"          "HORMAD1"         
## [124] "IL12RB2"          "TCAM1P"           "SLC15A1"         
## [127] "GABBR2"           "VGLL1"            "ROPN1"           
## [130] "AMY1A"            "EN1"              "MARCO"           
## [133] "LOC84740"         "A2ML1"            "ZIC1"            
## [136] "CA9"              "AQP5"             "FABP7"           
## [139] "MSLN"             "MUC16"            "ELF5"            
## [142] "PRAME"
```
While in the original dataset with only clinical data we had 817 patients, keeping only the patients that have overexpressed data for all that genes left us with a final amount of 270. 

Until the Cox Analysis we will only use the clinical data and put this aside by now.

-------

#3. Clasification with Hierarchical clustering {#hierarchical}

As an extra, using another learned method in  *Lab 1*, we will remove the variable living status (`cens`) from the dataframe and see if the clusters of Hierarchical clustering are able to predict if the patient is censored or dead basing on the rest of the covariates.

"Prediction" using  Hierarchical Clustering Methods is something that has already been done before as wee can see in some papers like [Dengue Prediction Using Hierarchical Clustering Methods](https://link.springer.com/chapter/10.1007/978-3-319-91800-6_11) or many others.


```r
make.clusters<-function(set,cluster.type, k) {

  d<-dist(set, method="euclidean") 
  clusters<-hclust(d, method=cluster.type)
  plot(clusters, labels=set$patient_id, hang = -1)
  rect.hclust(clusters,k,border=2:(2+k))
  set.cut<-cutree(clusters, k) 
  table(set.cut)
}

make.clusters(BC.no.cens,"complete",2)
```

![](paper_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```
## set.cut
##   1   2 
## 741  76
```

```r
table(BC$cens)
```

```
## 
##   0   1 
##  79 738
```

If we see the number of patients in each of the clusters is very similar to the number of patients for each level in the real status.


```r
hcut.cmp<-hcut(BC.no.cens, k=2, hc_method = "complete")
fviz_silhouette(hcut.cmp)
```

```
##   cluster size ave.sil.width
## 1       1  741          0.56
## 2       2   76          0.69
```

![](paper_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

The red color represent the patients which have been predicted to be censured, while the blue are the death ones. 
Excepting a small group, all the censured patiens have a "security" of more than 0.5 index.
The "death" cluster has a high silhouette index  (0.7). 

For verifying the results are as good as they seem, we will check the **accuracy** ($\frac{Coincidences}{Total}$) comparing the if the patient is in the cluster of censured or dead and comparing to the real value of the covariate `cens`.


```r
clusters<-hclust(dist(BC.no.cens, method="euclidean") , method="complete")
set.cut<-cutree(clusters, 2) -1 # 0 aqui son 1 de cens y 1 son 0 de cens
inverse.set.cut<-vector()
for(val in 1:length(set.cut)) {
  ifelse(set.cut[val]==0, inverse.set.cut[val]<-1, inverse.set.cut[val]<-0)
}
equals<-sum(inverse.set.cut==BC$cens)
different<-sum(inverse.set.cut!=BC$cens)
acc<-equals/(equals+different)
acc
```

```
## [1] 0.8372093
```

$0.83$ is a really good value, the experiment have worked really well, and at least in this case, Hierarchical clustering seems to be a good method to predict the status of survival patients basing of the rest of covariates.


### 3.1. SOM: Distribution of covariates in the clusters

Once we have the clusters, using Self Organization maps (also from *Lab 1*) we can see how the covariates are distributed in each cluster.


```r
numeric.BC<-BC.no.cens[-1]
factor.cols<-c("menostat", "ER", "PR", "HER2", "tgrade","Luminal", "NodalStatus", "NodeStatus", "NodalStatus2", "ClusteredAge")

for (col in factor.cols) numeric.BC[,col]<-as.numeric(numeric.BC[,col])

library(kohonen)
set.seed(100)
som_model <-som(scale(numeric.BC), grid=somgrid(7, 6, "hexagonal"), rlen=200, alpha =c(0.08,0.01))
#plot
som.hc = cutree(hclust(dist(som_model$codes[[1]])), 2)
plot(som_model, type="codes", bgcol=rainbow(2)[som.hc])
#cluster boundaries
add.cluster.boundaries(som_model, som.hc)
```

![](paper_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

---------



# 4. Cox Model and Paper comparison


## 4.1. Best covariates for Cox Model {#best_parameters_cox}

We will try to find which is the best combination of clinical covariates for Cox model For that we will get a previous idea making an exhaustive search and then we will make the final decision seeing how good the Cox model predictions are.

### 4.1.1. Exhaustive Search

Using the function of *Lab 3:Survival Analysis* we will try to find the best combination of covariates for our Cox Model with this specific data.




```r
library(FSelector)

exhaustiveSearch <-function(lista.vbles) {
  
  pvalue.acumulado<-list()
  
  for(nivel in 1:length(lista.vbles)) { #cada nivel (#vbles=#nivel)
    best.pvalue<-1
    for(v in lista.vbles) { #vbles de cada nivel
      ifelse(nivel==1, att<-v, att<-c(v, lista.acumulada))
      ec <-as.simple.formula(att,"Surv(time,cens)")
      sum<-summary(coxph(ec, data=BC))
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

es<-exhaustiveSearch(c("age", "menostat","tgrade","pnodes", "Luminal", "NodalStatus", "NodalStatus2", "NodeStatus", "ClusteredAge"))

df.order<-es[order(es$pvalue_acumulado),]
df.order
```

```
##   vble acumulada pvalue_acumulado
## 3         tgrade     5.956164e-09
## 4         pnodes     9.276238e-09
## 5            age     2.143009e-08
## 2       menostat     3.034877e-08
## 6   ClusteredAge     5.148877e-08
## 7     NodeStatus     1.249569e-07
## 8    NodalStatus     2.340871e-07
## 9        Luminal     1.984402e-06
## 1   NodalStatus2     6.423682e-06
```

The results are similar to the ones obtained in *Lab 3*, ponting out that there are some common important covariates from explaining survival in different types of cancer.

Now we have the model, before using it, we have to make sure it verifys the assumption of proportional risks: “If an individual has a risk of death at some initial time point that is twice as high as that of another individual, then at all later times the risk of death remains twice as high.” For that purpose, we’re going to use the function cox.zph(): the model will be accepted if $p−value≥0.05$.


```r
res.cox <- coxph(Surv(time,cens)~tgrade+pnodes+age, data=BC)
cox.zph(res.cox, transform="km", global=TRUE) 
```

```
##                rho    chisq       p
## tgradeII  -0.09997  7.36577 0.00665
## tgradeIII -0.05426  2.26663 0.13219
## tgradeIV  -0.05117  1.95267 0.16230
## pnodes    -0.00139  0.00129 0.97135
## age       -0.05796  2.29740 0.12959
## GLOBAL          NA 10.05718 0.07363
```

The global value is $0.07$ ≥$0.05$, therefore the model is accepted.



```r
summ<-summary(res.cox)
summ
```

```
## Call:
## coxph(formula = Surv(time, cens) ~ tgrade + pnodes + age, data = BC)
## 
##   n= 817, number of events= 738 
## 
##               coef exp(coef) se(coef)      z Pr(>|z|)    
## tgradeII   0.11306   1.11970  0.09985  1.132   0.2575    
## tgradeIII  0.06828   1.07066  0.14276  0.478   0.6325    
## tgradeIV  -0.97998   0.37532  0.51225 -1.913   0.0557 .  
## pnodes     0.02769   1.02808  0.01100  2.517   0.0118 *  
## age        0.01133   1.01140  0.00283  4.005 6.21e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##           exp(coef) exp(-coef) lower .95 upper .95
## tgradeII     1.1197     0.8931    0.9207     1.362
## tgradeIII    1.0707     0.9340    0.8093     1.416
## tgradeIV     0.3753     2.6644    0.1375     1.024
## pnodes       1.0281     0.9727    1.0061     1.050
## age          1.0114     0.9887    1.0058     1.017
## 
## Concordance= 0.577  (se = 0.012 )
## Rsquare= 0.037   (max possible= 1 )
## Likelihood ratio test= 30.67  on 5 df,   p=1e-05
## Wald test            = 30.81  on 5 df,   p=1e-05
## Score (logrank) test = 31.15  on 5 df,   p=9e-06
```

```r
summ$logtest[3]
```

```
##       pvalue 
## 1.090395e-05
```

```r
#plot(survfit(res.cox), xlab="Days", ylab="Portion not Rearrested")
```

The $p-value$ is really good, let's see how good it's predicting.

###  4.1.2. Prediction using Cox

Using the **three first covariates of the ranking**:


```r
predictions<-as.vector(sort(exp(-predict(res.cox, type=c("expected"))),decreasing=T))
survfit.values<-sort(survfit(formula=Surv(time,cens)~tgrade+pnodes+age, data=BC)$surv, decreasing=T)
all.times<-sort(BC$time)
```




```r
plot(all.times, predictions, ylim=c(0,1),col="red",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(all.times,survfit.values, ylim=c(0,1),col="blue",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("real","predicted"),lty=c(1,1),col=c(2,4))
```

![](paper_files/figure-html/unnamed-chunk-36-1.png)<!-- -->

The predicted survival values are quite bad, at the begining they follow the curve more or less, but at time 40 the probabilities start to be 0, while the real value is still at 0.5.

Let`s see if the predictions get better **using all the covariates**.



```r
all.cox <- coxph(Surv(time,cens)~., data=BC)
predictions<-as.vector(sort(exp(-predict(all.cox, type=c("expected"))),decreasing=T))
survfit.values<-sort(survfit(formula=Surv(time,cens)~1, data=BC)$surv, decreasing=T)
```




```r
plot(all.times, predictions, ylim=c(0,1),col="red",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(all.times,survfit.values, ylim=c(0,1),col="blue",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("real","predicted"),lty=c(1,1),col=c(2,4))
```

![](paper_files/figure-html/unnamed-chunk-39-1.png)<!-- -->

With all the covariates for Cox model the predictions are a little bit better than with only the three first ranked of exhaustive search but still not good.

Let's try an intermediate solution **adding the next covariates in the ranking**.



```r
cox2<-coxph(Surv(time,cens)~tgrade+pnodes+menostat+ClusteredAge+NodeStatus, data=BC)
predictions<-as.vector(sort(exp(-predict(cox2, type=c("expected"))),decreasing=T))
survfit.values<-sort(survfit(formula=Surv(time,cens)~tgrade+pnodes+menostat+ClusteredAge+NodeStatus, data=BC)$surv, decreasing=T)
```




```r
plot(all.times, predictions, ylim=c(0,1),col="red",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(all.times,survfit.values, ylim=c(0,1),col="blue",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("real","predicted"),lty=c(1,1),col=c(2,4))
```

![](paper_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

```r
sum<-summary(cox2)
p.value<-summ$logtest[3]
```

This is the best solution to the moment. We will keep the covariates **tgrade+pnodes+menostat+ClusteredAge+NodeStatus** for training the Cox model since they have demonstrate to have the best performance for this dataset.


## 4.2. Paper values using expression data {#comp_cox_paper}

We will try the Cox model with the selected best variables in the rpevious point and all the genes of the dataset made in [2.3.Expression Data](#expression_data) 





```r
#Formula with all gene names
gene.names<-names(BC.expression)[17:142]
gene.names.form<-paste(gene.names, collapse="+")
vbles.names.form<-paste(c("tgrade","pnodes", "menostat", "ClusteredAge", "NodeStatus"), collapse="+")
all.cov<-paste(c(vbles.names.form, gene.names.form), collapse="+")
forml=as.simple.formula(all.cov,"Surv(time,cens)")
all.times.expr<-sort(BC.expression$time)
#Cox using all genes
cox2.exp<-coxph(formula=forml, data=BC.expression)
predictions<-as.vector(sort(exp(-predict(cox2.exp, type=c("expected"))),decreasing=T))
survfit.values<-sort(survfit(formula=Surv(time,cens)~tgrade+pnodes+menostat+ClusteredAge+NodeStatus, data=BC.expression)$surv, decreasing=T)
```











```r
plot(all.times.expr, predictions, ylim=c(0,1),col="red",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
par(new=TRUE)
plot(all.times.expr ,survfit.values, ylim=c(0,1),col="blue",main="Survival curves using all the data",xlab="Time",ylab="Survival probability", type="l")
legend("topright",c("predicted using all genes","real"),lty=c(1,1),col=c(2,4))
```

![](paper_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

The Cox model using all genes covariates is the best one predicting till the moment.

#### Comparison 

About the parameters mentioned in [1.About the paper](#paper): *Cox Coefficient, p-value, FDR, Rank, Median Expresion* and *Mean Expresion*, we will not compare all of them. *FDR* (False discovery rate correction) and *RANK* is rank of the correlations for the different cancers. We can't calculate them cause we have only use one type of cancer.

The rest of parameters: 

- Cox Coefficient: `coef` colum, gives an idea of the impact of that gene in the rest of the model.

- Median Expression: Median expression of each gene in a specific cancer type.

- Mean Expression: Mean expression of each gene in a specific cancer type.

- P-Value.


We have the Cox Coefficients and the $p-values$ for all the genes in our Cox model:


```r
cox2.exp
```

```
## Call:
## coxph(formula = forml, data = BC.expression)
## 
##                               coef  exp(coef)   se(coef)      z        p
## tgradeII                 2.524e-01  1.287e+00  3.091e-01  0.817 0.414113
## tgradeIII                4.891e-01  1.631e+00  5.608e-01  0.872 0.383073
## tgradeIV                -2.579e+00  7.584e-02  1.448e+00 -1.781 0.074960
## pnodes                   8.290e-02  1.086e+00  3.917e-02  2.117 0.034287
## menostatPre              3.621e-01  1.436e+00  2.804e-01  1.291 0.196580
## ClusteredAge(60,90]      1.082e-01  1.114e+00  3.157e-01  0.343 0.731705
## NodeStatusNode positive  1.677e-01  1.183e+00  2.764e-01  0.607 0.544022
## AGR3                     8.934e-05  1.000e+00  9.756e-05  0.916 0.359827
## DNALI1                  -9.608e-04  9.990e-01  3.737e-04 -2.571 0.010132
## CAPN8                    2.520e-04  1.000e+00  4.104e-04  0.614 0.539135
## C6orf97                  2.523e-04  1.000e+00  4.939e-04  0.511 0.609446
## GPR77                   -9.693e-04  9.990e-01  2.501e-03 -0.388 0.698340
## SERPINA11                9.119e-05  1.000e+00  2.657e-04  0.343 0.731432
## FSIP1                    5.363e-04  1.001e+00  3.797e-04  1.412 0.157841
## TFF3                     2.508e-05  1.000e+00  2.920e-05  0.859 0.390410
## CACNA2D2                 5.467e-04  1.001e+00  4.943e-04  1.106 0.268764
## SYTL5                    9.638e-06  1.000e+00  5.924e-04  0.016 0.987019
## FLJ45983                -1.073e-03  9.989e-01  5.797e-04 -1.851 0.064139
## TFF1                    -4.869e-05  1.000e+00  2.660e-05 -1.831 0.067157
## PRR15                   -5.977e-04  9.994e-01  3.006e-04 -1.988 0.046757
## LOC145837               -1.382e-03  9.986e-01  1.148e-03 -1.204 0.228565
## DACH1                   -1.823e-04  9.998e-01  3.094e-04 -0.589 0.555713
## C1orf64                  4.743e-05  1.000e+00  3.077e-04  0.154 0.877503
## KCNK15                   2.865e-03  1.003e+00  9.773e-04  2.932 0.003371
## NAT1                     1.797e-04  1.000e+00  5.518e-05  3.256 0.001129
## ANKRD30B                 4.217e-04  1.000e+00  4.706e-04  0.896 0.370197
## MAPT                    -1.640e-04  9.998e-01  9.949e-05 -1.649 0.099188
## ERBB4                    3.588e-04  1.000e+00  1.681e-04  2.134 0.032805
## ESR1                     4.227e-05  1.000e+00  2.532e-05  1.669 0.095058
## SYT9                    -5.777e-04  9.994e-01  4.453e-04 -1.297 0.194526
## TSPAN1                   1.124e-04  1.000e+00  5.994e-05  1.875 0.060754
## NEK10                   -1.252e-03  9.987e-01  5.557e-04 -2.253 0.024250
## ZMYND10                  6.509e-04  1.001e+00  7.836e-04  0.831 0.406184
## DNAJC12                 -3.751e-05  1.000e+00  1.048e-04 -0.358 0.720486
## ABCC8                    1.496e-03  1.001e+00  5.675e-04  2.637 0.008372
## BCAS1                   -5.098e-04  9.995e-01  1.900e-04 -2.682 0.007310
## SLC44A4                  1.368e-05  1.000e+00  1.148e-04  0.119 0.905142
## CT62                    -2.224e-03  9.978e-01  1.721e-03 -1.292 0.196296
## AGR2                     2.132e-05  1.000e+00  1.216e-05  1.753 0.079657
## GP2                     -8.193e-05  9.999e-01  5.238e-05 -1.564 0.117775
## CHAD                     2.673e-05  1.000e+00  6.566e-05  0.407 0.683869
## PGR                     -1.674e-05  1.000e+00  1.640e-05 -1.021 0.307387
## GRIK3                   -2.297e-04  9.998e-01  2.741e-04 -0.838 0.402108
## DEGS2                   -3.959e-05  1.000e+00  1.489e-04 -0.266 0.790310
## SCUBE2                  -1.433e-05  1.000e+00  1.370e-05 -1.046 0.295457
## LOC100130148             6.232e-03  1.006e+00  6.838e-03  0.911 0.362090
## AKR7A3                   8.499e-05  1.000e+00  1.758e-04  0.483 0.628810
## TTC36                   -2.328e-03  9.977e-01  2.148e-03 -1.084 0.278548
## PTPRT                    1.564e-04  1.000e+00  9.019e-05  1.734 0.082972
## ADAMTS15                -7.648e-06  1.000e+00  4.057e-04 -0.019 0.984958
## GATA3                    4.136e-05  1.000e+00  3.491e-05  1.185 0.236181
## IGFALS                  -3.407e-04  9.997e-01  9.934e-04 -0.343 0.731644
## GREB1                    1.278e-04  1.000e+00  7.033e-05  1.818 0.069131
## ANKRD30A                -1.995e-04  9.998e-01  1.416e-04 -1.409 0.158761
## CA12                     7.067e-06  1.000e+00  1.402e-05  0.504 0.614115
## FOXA1                   -4.256e-06  1.000e+00  8.905e-05 -0.048 0.961883
## MLPH                     5.322e-05  1.000e+00  7.423e-05  0.717 0.473419
## TMPRSS6                 -1.857e-04  9.998e-01  9.166e-04 -0.203 0.839487
## AFF3                    -1.042e-04  9.999e-01  9.031e-05 -1.154 0.248376
## F7                       4.802e-04  1.000e+00  8.390e-04  0.572 0.567093
## LOC100128977            -2.096e-03  9.979e-01  7.403e-03 -0.283 0.777045
## C20orf114                2.890e-05  1.000e+00  2.835e-05  1.019 0.308107
## PGLYRP2                 -1.171e-03  9.988e-01  5.503e-04 -2.128 0.033329
## TPRG1                   -1.962e-04  9.998e-01  1.103e-04 -1.779 0.075253
## RGS22                    1.488e-03  1.001e+00  5.991e-04  2.483 0.013032
## GRPR                    -1.723e-04  9.998e-01  1.699e-04 -1.014 0.310504
## AR                       3.811e-04  1.000e+00  3.325e-04  1.146 0.251728
## SPDEF                   -2.989e-05  1.000e+00  9.781e-05 -0.306 0.759944
## CYP2B7P1                -8.181e-06  1.000e+00  2.543e-05 -0.322 0.747663
## SERPINA5                 1.053e-05  1.000e+00  7.916e-05  0.133 0.894168
## LOC728606               -6.295e-04  9.994e-01  1.059e-03 -0.594 0.552183
## GFRA1                   -1.037e-05  1.000e+00  2.067e-05 -0.502 0.615963
## EEF1A2                  -1.500e-05  1.000e+00  5.954e-05 -0.252 0.801085
## SLC7A2                  -4.900e-05  1.000e+00  2.921e-05 -1.677 0.093484
## TUBA3E                  -9.981e-04  9.990e-01  1.474e-03 -0.677 0.498289
## LRG1                    -1.174e-04  9.999e-01  1.219e-04 -0.963 0.335465
## CYP4B1                   3.043e-05  1.000e+00  7.569e-05  0.402 0.687656
## CYP4Z2P                 -4.706e-04  9.995e-01  2.621e-04 -1.795 0.072599
## KCNC2                   -3.658e-04  9.996e-01  2.272e-04 -1.610 0.107339
## MS4A8B                   3.102e-03  1.003e+00  8.610e-04  3.603 0.000315
## SYT1                    -1.966e-04  9.998e-01  8.336e-05 -2.359 0.018347
## EMX1                     1.115e-03  1.001e+00  2.388e-03  0.467 0.640448
## HMGCS2                  -1.102e-05  1.000e+00  4.511e-05 -0.244 0.806966
## C1orf173                -9.036e-04  9.991e-01  9.201e-04 -0.982 0.326095
## LOC389033               -2.152e-03  9.979e-01  2.688e-03 -0.801 0.423316
## SERPINA6                 5.593e-05  1.000e+00  4.151e-05  1.347 0.177887
## TMC5                    -9.873e-06  1.000e+00  9.946e-05 -0.099 0.920920
## HEPACAM2                -6.873e-04  9.993e-01  7.835e-04 -0.877 0.380395
## RIMS4                    7.343e-05  1.000e+00  1.810e-04  0.406 0.685001
## PHGR1                   -2.008e-03  9.980e-01  1.217e-03 -1.650 0.098863
## CHST8                    3.045e-04  1.000e+00  3.827e-04  0.796 0.426272
## CEACAM6                 -1.792e-06  1.000e+00  2.308e-05 -0.078 0.938120
## CEACAM5                  9.817e-05  1.000e+00  1.435e-04  0.684 0.494073
## ABCC11                  -3.994e-05  1.000e+00  5.681e-05 -0.703 0.482060
## SLC30A8                 -1.222e-05  1.000e+00  2.109e-05 -0.580 0.562228
## BPIL1                   -3.202e-04  9.997e-01  5.864e-04 -0.546 0.584979
## CPB1                     2.875e-06  1.000e+00  1.489e-06  1.930 0.053554
## KCNJ3                    6.141e-05  1.000e+00  2.093e-04  0.293 0.769167
## C10orf82                -2.573e-03  9.974e-01  1.166e-03 -2.206 0.027370
## CDC20B                   4.223e-04  1.000e+00  3.419e-04  1.235 0.216771
## SCGB2A2                  1.062e-06  1.000e+00  3.701e-06  0.287 0.774099
## DHRS2                    2.107e-07  1.000e+00  1.285e-05  0.016 0.986925
## ABCC13                   1.052e-04  1.000e+00  1.054e-03  0.100 0.920524
## BMPR1B                  -5.117e-06  1.000e+00  2.978e-05 -0.172 0.863548
## PIP                      7.643e-06  1.000e+00  6.934e-06  1.102 0.270326
## SCGB1D2                  1.182e-06  1.000e+00  1.355e-05  0.087 0.930485
## VSTM2A                  -2.373e-04  9.998e-01  1.184e-04 -2.004 0.045080
## CST9                     4.507e-04  1.000e+00  2.816e-04  1.600 0.109539
## CLEC3A                   5.571e-06  1.000e+00  1.254e-05  0.444 0.656992
## FZD9                     9.407e-04  1.001e+00  6.507e-04  1.446 0.148313
## ART3                    -3.533e-04  9.996e-01  2.695e-04 -1.311 0.189952
## CXorf49B                -1.981e-03  9.980e-01  8.537e-04 -2.320 0.020322
## CXorf61                 -1.797e-04  9.998e-01  1.127e-03 -0.159 0.873297
## HORMAD1                  3.792e-04  1.000e+00  4.657e-04  0.814 0.415488
## IL12RB2                  9.222e-04  1.001e+00  2.766e-04  3.334 0.000857
## TCAM1P                  -1.061e-03  9.989e-01  1.379e-03 -0.769 0.441744
## SLC15A1                  3.929e-04  1.000e+00  2.240e-04  1.754 0.079388
## GABBR2                  -3.321e-04  9.997e-01  3.980e-04 -0.834 0.404039
## VGLL1                   -2.722e-06  1.000e+00  7.620e-05 -0.036 0.971502
## ROPN1                   -5.077e-05  9.999e-01  2.160e-04 -0.235 0.814199
## AMY1A                    3.235e-04  1.000e+00  1.844e-04  1.754 0.079407
## EN1                     -8.127e-05  9.999e-01  1.544e-04 -0.526 0.598589
## MARCO                   -9.700e-04  9.990e-01  4.559e-04 -2.127 0.033379
## LOC84740                 2.704e-04  1.000e+00  1.923e-04  1.406 0.159647
## A2ML1                    7.111e-05  1.000e+00  6.346e-05  1.120 0.262539
## ZIC1                     1.364e-03  1.001e+00  6.047e-04  2.255 0.024133
## CA9                     -1.260e-04  9.999e-01  1.658e-04 -0.760 0.447455
## AQP5                     1.126e-04  1.000e+00  7.276e-05  1.547 0.121853
## FABP7                    9.810e-06  1.000e+00  7.186e-06  1.365 0.172191
## MSLN                    -8.560e-05  9.999e-01  4.280e-05 -2.000 0.045523
## MUC16                   -3.781e-05  1.000e+00  5.958e-05 -0.635 0.525726
## ELF5                    -1.149e-04  9.999e-01  3.753e-05 -3.061 0.002203
## PRAME                    5.829e-05  1.000e+00  9.268e-05  0.629 0.529377
## 
## Likelihood ratio test=222  on 133 df, p=2.035e-06
## n= 270, number of events= 242
```


For comparing the results with the paper ones, we have to choose a gene, for example *"MARCO"*.

[MARCO paper analysis](http://www.oncolnc.org/search_results/?q=MARCO+)



```r
mean(BC.expression$MARCO)
```

```
## [1] 134.9812
```

```r
median(BC.expression$MARCO)
```

```
## [1] 20.96855
```

|     | Cox Coefficient | P-value | Median Expression | Mean Expression |
|------|------------------|--------|--------------------|-----------------|
|oncolnc| -1.6e-03      | 8.40e-01| 14.11             | 110.97           |
|ours| -9.700e-04      |3.3379e-01|  20.9685        | 134.98 |

We will compare some others:

- **CST9**: Is not avaliable in oncolnc.

- **KCNJ3**: Only avaliable for 8 cancers. (In the page it's said: *Missing cancers? Some genes don't meet the expression cutoff for the analysis, or there was no expression data for that cancer for that gene.*)
 
Once saw that not all of our genes are avaliable in the application, we will search for ones with a higher Cox Coefficient.

- **GPR77** they have it's synonim *C5AR2*. 


```r
mean(BC.expression$GPR77)
```

```
## [1] 97.07879
```

```r
median(BC.expression$GPR77)
```

```
## [1] 83.7962
```

|     | Cox Coefficient | P-value | Median Expression | Mean Expression |
|------|------------------|--------|--------------------|-----------------|
|oncolnc| -3.2e-03	     | 8.40e-01| 95.06            |114.57          |
|ours| -9.693e-04     |6.98340e-01| 83.79       | 97.07|

- **EN1**


```r
mean(BC.expression$EN1)
```

```
## [1] 450.3667
```

```r
median(BC.expression$EN1)
```

```
## [1] 39.0473
```

|     | Cox Coefficient | P-value | Median Expression | Mean Expression |
|------|------------------|--------|--------------------|-----------------|
|oncolnc| -8.127e-05     | 5.998589e-01| 34.19           |301.04        |
|ours| -9.693e-04     |6.98340e-01|  39.0473      |  450.3667|

The results are surprisingly very similar. The higher (but not a lot) difference is found, as expected, in the comparision oc Coex Coefficients. This is because this value depends on the rest of the dataset, and they have used a different one. 

We can see some important differences of *oncolnc* and our dataset with a summary of their data they provided (in [Table 1](https://peerj.com/articles/cs-67/)):


```r
nrow(BC.expression) #number of patients
```

```
## [1] 270
```

```r
mean(BC.expression$age) #mean age at diagnosis
```

```
## [1] 57.88148
```

```r
sum(BC.expression$cens==1) #events
```

```
## [1] 242
```

```r
length(gene.names) #number of genes
```

```
## [1] 126
```

|     | Cancer | Patients | Male/Female |Age at diagnosis | Events | Genes |
|------|------------------|-------------|------------------|------|--------|
|oncolnc| BRCA| 1.006     | 11/995      |  58.34             |135 | 16,602 |
|ours| BRCA| 270     | *--unkonown--*  |  57.88            |242 | 126|



The datasets seem to be quite similar, excepting we have much more events and they have much more patients.

That can explain the sightly differences, but in general, we have obtain good similarities.



