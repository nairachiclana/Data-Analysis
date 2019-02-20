#!/bin/bash

#color options
red=`tput setaf 1`
blue=`tput setaf 6`
green=`tput setaf 2`
reset=`tput sgr0`

install_packages=false
volcano_plot=false
survival_curves=false

while getopts "ivl:p:e:x:sf:c:r:k:t:a:b:" opt; do
    case "$opt" in

    i)
		install_packages=true
		;;
    v)
        volcano_plot=true
        ;;

    l)
		lcf_volcano=${OPTARG}
		;;
	p)
		pval_volcano=${OPTARG}
		;;
	e)
		lcf_sobreexpresed_data=${OPTARG}
		;;
	x)
		pval_sobreexpresed_data=${OPTARG}
		;;

	s)
		survival_curves=true
		;;
	f)
		subset_filter=${OPTARG}
		;;
	c)
		subset_cox=${OPTARG}
		;;
	r)
		lr_k=${OPTARG}
		;;
	k)
		cv_k=${OPTARG}
		;;
	t)
		type_expression==${OPTARG}
		;;
	a)
		lfc_fa==${OPTARG}
		;;
	b)
		pval_fa==${OPTARG}
		;;

    
    
    esac
done

shift $((OPTIND-1))

#Install de packegaes
if  $install_packages 
then
	RScript 0_install_packages.R 
	echo "${green} Se han instalado los paquetes ${reset}"
fi

#1.DATA

#1.1.Clinical data
echo  "${blue} Matrix of clinical data ${reset}"
RScript 1_1_Clinical_Data.R 

#1.2.Expression data
echo  "${blue} Matrix of expression data ${reset}"
RScript 1_2_Expression_Data.R 

#VolcanoPlot
if  $volcano_plot  
then
	echo "${blue} Procceding to make the volcano plot ${reset}": 
	RScript VolcanoPlot.R $lcf_volcano $pval_volcano
	echo "${green} Se ha guardado el volcano plot  en la carpeta de origen ${reset}"
fi


#1.3.Join Sobreexpressed and clinical data
echo "${blue} Procceding to join clinical and expression data ${reset}": 
RScript 1_3_Join_clinical_sobreexpresed.R $lcf_sobreexpresed_data $pval_sobreexpresed_data
echo "${green} Se ha realizado el join ${reset}"

#2.DATA 
#2.1.Stratification with clustering 
echo "${blue} Procceding to stratify the covariates ${reset}": 
RScript 2_1_Data_stratification.R 

#2.1.Survival with ggplot 
if  $survival_curves 
then
	echo "${blue} Procceding to do survival curves ${reset}": 
	RScript 2_2_Survival_Curves.R 
	echo "${green} Se han guardado los plots en la carpeta results ${reset}"
fi

echo "${blue} Procceding to do variable selection ${reset}"
RScript 2_3_Cox_Analysis.R  $subset_filter $subset_cox


#3.PREDICTIVE MODELS

#3.1.Data replication
echo "${blue} Procceding to replicate the data ${reset}": 
RScript 3_1_Replicate_Data.R  
echo "${green} Se ha replicado el dataset ${reset}"

#3.2. Lasso and Ridge Regression
echo "${blue} Procceding to execute Lasso and Ridge ${reset}":
RScript 3_2_Lasso_Ridge.R  $lr_k

#3.3.Predictors
echo "${blue} Procceding to execute GLM, DT, NNET ${reset}":
RScript 3_3_GLM_DT_NNET.R  $cv_k
echo "${green} Se ha guardado el Boxplot en la carpeta results ${reset}"

#4.FUNCTIONAL ANALYSIS
echo "${blue} Procceding to execute functional analysis ${reset}":
RScript 4_FunctionalAnalysis.R  $type_expression $lfc_fa $pval_fa
echo "${green} Se ha guardado el dotplot en la carpeta origen ${reset}"



