### validation_code.R

### Author: Valeria Vitelli
### Contributers: Andrew Reiner, Alvaro Köhn-Luque

### Validations of the models in Xie et al., 2020; Zhang et al., 2020; Allenbach et al., 2020
### See validation_README.txt in this directory for the format of the
### input file 'data.xlsx'

### rm(list=objects())

library(readxl)
library(dplyr)
library(ROCit)
library(VIM)
source('./validation_functions.R')

ValidateOnCohort <- FALSE # should real data be used? (alternative: simulate)
getSummaries <- TRUE      # should final results summaries and plots be produced?

## LOAD DATA

if(ValidateOnCohort) data <- read_excel('./data.xlsx') else data <- simulate.data(200)$data
head(data)

# Confirm that validation dataset contains required features/outcomes 
required <- simulate.data()$features
if(sum(is.na(match(required, names(data))))>0) 
  print('Error: not all variables are provided in the dataset; please check.')

#############################################################
##### Xie model validation:
#############################################################

## prepare the data
variable_Xie <- c("death","age","lymp","ld","spo2")
data_Xie <- data[,variable_Xie]

## if needed: data imputation, default = k-NN (new data are called "data_imp")
if(sum(is.na(data_Xie))>0) data_imp <- impute.mydata(data_Xie) else data_imp <- data_Xie

## Xie et al. model computations
resXie <- validate.Xie(data_imp)

## recalibration, default = adjustment of the intercept (baseline risk/hazard)
## (done for adjusting the outcome frequency across cohorts)
corrXie <- recalibrate(resXie)
resXie <- validate.Xie(data_imp, recalib = TRUE, correction = corrXie)


#############################################################
##### Zhang model validation:
#############################################################

## prepare the data
variable_Zhang <- c("death","poor_outcome","age","gender","crp",
                "neut","lymp","plt","crea")
data_Zhang <- data[,variable_Zhang]

## if needed: data imputation, default = k-NN (new data are called "data_imp")
if(sum(is.na(data_Zhang))>0) data_imp <- impute.mydata(data_Zhang) else data_imp <- data_Zhang

## Zhang et al. model computations
resZhang <- validate.Zhang(data_imp)

## recalibration, default = adjustment of the intercept (baseline risk/hazard)
## (done for adjusting the outcome frequency across cohorts)
if(dim(resZhang$out)[2]>1){
  corrZhang <- c(recalibrate(list(out=resZhang$out[,1],pred=resZhang$pred[,1])), # outcome death
                 recalibrate(list(out=resZhang$out[,2],pred=resZhang$pred[,2])) ) # poor outcome
}else corrZhang <- recalibrate(resZhang)
resZhang <- validate.Zhang(data_imp, recalib = TRUE, correction = corrZhang)


#############################################################
##### Allenbach model validation:
#############################################################

## prepare the data
variable_Allenbach <- c("deathICU","age","WHOscale","lymp","crp")
data_Allenbach <- data[,variable_Allenbach]

## if needed: data imputation, default = k-NN (new data are called "data_imp")
if(sum(is.na(data_Allenbach))>0) data_imp <- impute.mydata(data_Allenbach) else data_imp <- data_Allenbach

## Allenbach et al. model computations
resAllenbach <- validate.Allenbach(data_imp)

## recalibration, default = adjustment of the intercept (baseline risk/hazard)
## (done for adjusting the outcome frequency across cohorts)
corrAllenbach <- recalibrate(resAllenbach)
resAllenbach <- validate.Allenbach(data_imp, recalib = TRUE, correction = corrAllenbach)

#############################################################

## SAVE FINAL RESULTS
save(resXie, resZhang, resAllenbach, file = './validation_results.Rdata')

## if needed, SUMMARIZE & VISUALIZE final results
if(getSummaries){
  
  pdf('./validation_ROC.pdf', height = 10, width = 10)
  par(mfrow = c(2,2))
  
  print('Results of validation of Xie et al., 2020:')
  rocXie <- rocit(resXie$pred,resXie$out)
  summary(rocXie)
  plot(rocXie, main = 'Xie et al., 2020')
  
  print('Results of validation of Allenbach et al., 2020:')
  rocAllenbach <- rocit(resAllenbach$pred,resAllenbach$out)
  summary(rocAllenbach)
  plot(rocAllenbach, main = 'Allenbach et al., 2020')
  
  print('Results of validation of Zhang et al., 2020:')
  if(dim(resZhang$out)[2]>1){
    print('Outcome death:')
    rocZhangD <- rocit(resZhang$pred[,1],resZhang$out[,1])
    summary(rocZhangD)
    plot(rocZhangD, main = 'Zhang et al., 2020 - outcome death')
    print('Poor outcome:')
    rocZhangP <- rocit(resZhang$pred[,2],resZhang$out[,2])
    summary(rocZhangP)
    plot(rocZhangP, main = 'Zhang et al., 2020 - poor outcome')
  }else{
    rocZhang <- rocit(as.vector(resZhang$pred),as.vector(resZhang$out))
    summary(rocZhang)
    plot(rocZhang, main = 'Zhang et al., 2020')
  }
  
  dev.off()
}

