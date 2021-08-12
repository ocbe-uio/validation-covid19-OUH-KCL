### validation_functions.R

### Author: Valeria Vitelli
### Contributers: Andrew Reiner, Alvaro Köhn-Luque

### Support functions for validation_code.R

##-----------------------------
## function for data generation
##-----------------------------
## (if the user wishes *not* to validate a cohort)
simulate.data <- function(n=500){
  # outcomes simulated from Bernoulli with reasonable p
  outcomes <- cbind(rbinom(n,1,.1),rbinom(n,1,.5),rbinom(n,1,.2))
  # age ~ N(50,10^2), gender ~ Be(.5), all other features are ~N(0,1) 
  features <- cbind(round(rnorm(n,50,10),0),rbinom(n,1,.5),
                    matrix(rnorm(n*8),n,8))
  tmp <- cbind(outcomes, features)
  data <- data.frame(tmp)
  myvars <- c('death', 'poor_outcome', 'deathICU', 'age', 'gender', 'lymp', 
              'ld', 'spo2', 'crp', 'neut', 'plt','crea','WHOscale')
  names(data) <- myvars
  return(list(data=data, features=myvars))
}

##-----------------------------
## function for data imputation
##-----------------------------
impute.mydata <- function(data, method='kNN'){
  data <- kNN(data)
  return(data)
}

##--------------------------------------
## function for the model re-calibration
##--------------------------------------
recalibrate <- function(model, method='method1'){
  # 'model': a model estimated by validate.Xie, validate.Zhang, or validate.Allenbach
  
  obs.out <- mean(model$out)
  mean.pred <- mean(model$pred)
  correction.factor <- log((obs.out/(1-obs.out))/(mean.pred/(1-mean.pred)))
  return(correction.factor)
}

##------------------------------------------------
## function for the validation of Xie et al., 2020
##------------------------------------------------
validate.Xie <- function(data, recalib = FALSE, correction = NA){
  # 'recalib': boolean; does the model needs recalibration?
  # 'correction': scalar; ignored if recalib is FALSE, should be non-NA otherwise
  
  # check coherence
  if(recalib & is.na(correction))
    print('Error: if the model needs recalibration, correction factor should be provided; please check.')

  # store the complete outcomes/covariates explicitly, and apply data transformations
  outcome <- data_imp$death
  age <- data_imp$age
  LYM <- data_imp$lymp
  LYM[which(LYM<=0)] <- 1 # should not happen
  logLYM <- log(LYM) # this is actually used in the model
  LDH <- data$ld
  SPO <- data$spo2
  
  # model coefficients 
  # (Table 3, pag. 19 of Xie et al., "Optimism adjusted beta coefficients")
  age.coef <- 0.047
  LDH.coef <- 0.003
  LYM.coef <- -1.094
  SPO.coef <- -0.098
  intercept <- 4.559
  if(recalib)intercept <- intercept + correction
  
  # model fit
  lin.exp <- intercept + age.coef*age + LDH.coef*LDH + LYM.coef*logLYM + SPO.coef*SPO
  predicted.values <- as.vector(1/(1+exp(-lin.exp)))
  
  return(list(out = outcome, pred = predicted.values))
  
}


##--------------------------------------------------
## function for the validation of Zhang et al., 2020
##--------------------------------------------------
validate.Zhang <- function(data, which.outcome = 'both', recalib = FALSE, correction = NA){
  # 'which.outcome': character string; specifies which outcomes need to be validated.
  #                  Must be one of 'death', 'poor' or 'both'. Defaults to both outcomes.
  # 'recalib': boolean; does the model needs recalibration?
  #            (can be a 2-dim vector, in which case its elements refer to outcome 'death' and 'poor', respectively)
  # 'correction': scalar; ignored if recalib is FALSE, should be non-NA otherwise
  #             (has to be a 2-dim vector if recalibration is needed for both outcomes)
  
  # check coherence
  if(any(recalib) & all(is.na(correction)))
    print('Error: if the model needs recalibration, correction factor should be provided; please check.')
  if((which.outcome == 'both') & all(recalib) & (length(correction)==1))
    print('Error: if recalibration is needed for both outcomes, 
          correction factor should be specified for each model outcome; please check.')
  
  # store the complete outcomes/covariates explicitly, and apply data transformations
  outcomeDeath <- data_imp$death
  outcomePoor <- data_imp$poor_outcome
  age <- data_imp$age
  gender <- data_imp$gender
  crp <- data_imp$crp
  neut <- data_imp$neut
  lymp <- data_imp$lymp
  plt <- data_imp$plt
  crea <- data_imp$crea
  outcome <- predicted.values <- NULL
  
  if(which.outcome != 'poor'){# runs for 'death' and 'both'
    
    # model coefficients - death
    # (Supplementary Table 2, pag.10-11, Supplementary Material of Zhang et al.)
    age.coef <- 0.0142
    gender.coef <- -0.8035
    crp.coef <- 0.0138
    neut.coef <- 0.2249
    lymp.coef <- -1.3255
    plt.coef <- -0.0052
    crea.coef <- 0.0015
    intercept <- -3.5172
    if(recalib[1])intercept <- intercept + correction[1]
    
    # model fit
    lin.exp <- intercept + age.coef*age + gender.coef*gender + crp.coef*crp 
               + neut.coef*neut + lymp.coef*lymp + plt.coef*plt + crea.coef*crea
    predicted.values.death <- as.vector(1/(1+exp(-lin.exp)))
    
    # final set-up
    outcome <- cbind(outcome, outcomeDeath)
    predicted.values <- cbind(predicted.values, predicted.values.death)
    
  }

  if(which.outcome != 'death'){# runs for 'poor' and 'both'
    
    # model coefficients - poor outcome
    # (Supplementary Table 2, pag.10-11, Supplementary Material of Zhang et al.)
    age.coef <- 0.0281
    gender.coef <- -0.4723
    crp.coef <- 0.0082
    neut.coef <- 0.2575
    lymp.coef <- -1.1842
    plt.coef <- -0.0027
    crea.coef <- 0.0099
    intercept <- -4.1865
    if(recalib[length(recalib)])intercept <- intercept + correction[length(correction)]
    
    # model fit
    lin.exp <- intercept + age.coef*age + gender.coef*gender + crp.coef*crp 
               + neut.coef*neut + lymp.coef*lymp + plt.coef*plt + crea.coef*crea
    predicted.values.poor <- as.vector(1/(1+exp(-lin.exp)))

    # final set-up
    outcome <- cbind(outcome, outcomePoor)
    predicted.values <- cbind(predicted.values, predicted.values.poor)
    
  }
  

  return(list(out = outcome, pred = predicted.values))
  
}


##------------------------------------------------------
## function for the validation of Allenbach et al., 2020
##------------------------------------------------------
validate.Allenbach <- function(data, recalib = FALSE, correction = NA){
  # 'recalib': boolean; does the model needs recalibration?
  # 'correction': scalar; ignored if recalib is FALSE, should be non-NA otherwise
  
  # check coherence
  if(recalib & is.na(correction))
    print('Error: if the model needs recalibration, correction factor should be provided; please check.')
  
  # store the complete outcomes/covariates explicitly, and apply data transformations
  outcome <- data_imp$deathICU
  age <- data_imp$age
  age01 <- as.numeric(age > 60) # this is actually used in the model
  WHOscale <- data_imp$WHOscale
  lymp <- data_imp$lymp
  crp100 <- data_imp$crp
  crp <- crp100/100 # this is actually used in the model (CRP in 100mg/L)
  
  # model coefficients 
  # (Table S2, pag. 2, Supplementary Material of Allenbach et al.)
  age.coef <- 0.96
  WHO.coef <- 1.395
  lymp.coef <- -1.035
  crp.coef <- 0.49
  intercept <- -6.584
  if(recalib)intercept <- intercept + correction
  
  # model fit
  lin.exp <- intercept + age.coef*age01 + WHO.coef*WHOscale + lymp.coef*lymp + crp.coef*crp
  predicted.values <- as.vector(1/(1+exp(-lin.exp)))
  
  return(list(out = outcome, pred = predicted.values))
  
}
