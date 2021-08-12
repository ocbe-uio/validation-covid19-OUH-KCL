### validation_all_models.R
###
### Author: Valeria Vitelli
### Contributers: Andrew Reiner, Alvaro Köhn-Luque
###
### Validation of all models on the Oslo University Hospital Data
### For processing of other data, it is recoomended to use
### validation_code_London/validation_code.R
### 
###
### This script expects imputed data to be in subdirectory 'imputation'
###
### This script writes output to subdirectories 'Allenbach', 'Xie', 'Zhang', 'validation_plots_all_models', and expects them to exist



### rm(list=objects())

imputation_dir = 'imputation'

library(readxl)
library(dplyr)
library(ROCit)

####-----------------------------------------------
#### validation of all models on OSLO data:
#### prepare the data for the plotting / validation
####-----------------------------------------------

## Allenbach
##----------

which.impu.type <- c('singleKNN')

for(type in which.impu.type){
  
  switch(match(type, which.impu.type),
         {B <- 0; impu.type <- 'KNN'},
         {B <- 0; impu.type <- 'rf'},
         {B <- 9; impu.type <- 'br'},
         {B <- 9; impu.type <- 'gp'})
  
  predval.all <- NULL
  for(b in 0:B){
    
    data <- read_excel(paste(imputation_dir, '/Allenbach/validation_Allenbach_fill_',impu.type,'_',b,'.xlsx',sep=''))
    
    ## outcome
    outcome <- pull(data[,match('ICU or death by 14 days:', names(data))])    
    ## covariates
    age <- pull(data[,match('Alder_Ny', names(data))])
    age01 <- as.numeric(age > 60)

    WHOscale <- pull(data[,match('WHO scale', names(data))])
    LYM <- pull(data[,match('Lymfocytter', names(data))])
    CRP100 <- pull(data[,match('CRP', names(data))])
    CRP <- CRP100/100

    ## model coefficients (from Allenbach et al., Table S2)
    age.coef <- 0.96
    WHO.coef <- 1.395
    LYM.coef <- -1.035
    CRP.coef <- 0.49
    intercept <- -6.584
    
    ## predictions
    lin.exp <- intercept + age.coef*age01 + WHO.coef*WHOscale + LYM.coef*LYM + CRP.coef*CRP
    predicted.values <- as.vector(1/(1+exp(-lin.exp)))
    
    ## recalibration (method 1)
    obs.out <- mean(outcome)
    mean.pred <- mean(predicted.values)
    correction.factor <- log((obs.out/(1-obs.out))/(mean.pred/(1-mean.pred)))
    intercept <- intercept + correction.factor
    lin.exp <- intercept + age.coef*age01 + WHO.coef*WHOscale + LYM.coef*LYM + CRP.coef*CRP
    predicted.values <- as.vector(1/(1+exp(-lin.exp)))
    predval.all <- cbind(predval.all, predicted.values)

  }## end loop on imputation tables
  
  # save
  save(predval.all, outcome, file=paste('./Allenbach/validation_Allenbach_',type,'.RData',sep=''))
  
}## end loop on imputation method


## Xie
##----

which.impu.type <- c('singleKNN')

for(type in which.impu.type){
  
  switch(match(type, which.impu.type),
         {B <- 0; impu.type <- 'KNN'},
         {B <- 0; impu.type <- 'rf'},
         {B <- 9; impu.type <- 'br'},
         {B <- 9; impu.type <- 'gp'})
  
  predval.all <- NULL
  for(b in 0:B){
    
    ## read data
    data <- read_excel(paste(imputation_dir, '/Xie/validation_Xie_fill_',impu.type,'_',b,'.xlsx',sep=''))
    names(data) <- c('ID', 'death', 'age', 'lymp', 'log10lymph', 'spo2', 'ld')
    ## outcome
    outcome <- data$death
    
    ## covariates:
    # age
    age <- data$age
    # LDH
    LDH <- data$ld
    # Lymphocite
    LYM <- data$lymp
    LYM[which(LYM<=0)] <- 1
    logLYM <- log(LYM)
    # SPO2
    SPO <- data$spo2
    
    ## model coefficients (from Xie et al.)
    age.coef <- 0.047
    LDH.coef <- 0.003
    LYM.coef <- -1.094
    SPO.coef <- -0.098
    intercept <- 4.559
    
    ## predictions
    lin.exp <- intercept + age.coef*age + LDH.coef*LDH + LYM.coef*logLYM + SPO.coef*SPO
    predicted.values <- as.vector(1/(1+exp(-lin.exp)))
    
    ## recalibration (method 1)
    obs.out <- mean(outcome)
    mean.pred <- mean(predicted.values)
    correction.factor <- log((obs.out/(1-obs.out))/(mean.pred/(1-mean.pred)))
    intercept <- intercept + correction.factor
    lin.exp <- intercept + age.coef*age + LDH.coef*LDH + LYM.coef*logLYM + SPO.coef*SPO
    predicted.values <- as.vector(1/(1+exp(-lin.exp)))
    
    predval.all <- cbind(predval.all, predicted.values)
    
  }## end loop on imputation tables
  
  # save
  save(predval.all, outcome, file=paste('./Xie/validation_Xie_',type,'.RData',sep=''))
  
}## end loop on imputation method

## Zhang
##------

which.impu.type <- c('singleKNN')

for(type in which.impu.type){
  
  switch(match(type, which.impu.type),
         {B <- 0; impu.type <- 'KNN'},
         {B <- 0; impu.type <- 'rf'},
         {B <- 9; impu.type <- 'br'},
         {B <- 9; impu.type <- 'gp'})
  
  predval.all.death <- predval.all.poor <- NULL
  for(b in 0:B){
    
    ## read data - death outcome
    data <- read_excel(paste(imputation_dir, '/Zhang/validation_Zhang1_fill_',impu.type,'_',b,'.xlsx',sep=''))
    
    ## outcome
    outcomeDeath <- pull(data[,match('Death (inhospital)', names(data))])
    
    ## covariates
    age <- pull(data[,match('Alder_Ny', names(data))])
    gender <- pull(data[,match('sex', names(data))])
    crp <- pull(data[,match('CRP', names(data))])
    neut <- pull(data[,match('neutrophils', names(data))])
    lymp <- pull(data[,match('Lymfocytter', names(data))])
    plt <- pull(data[,match('platelets', names(data))])
    crea <- pull(data[,match('kreatinin', names(data))])
    
    ## model coefficients
    # (Supplementary Table 2, pag.10-11, Supplementary Material of Zhang et al.)
    age.coef <- 0.0142
    gender.coef <- -0.8035
    crp.coef <- 0.0138
    neut.coef <- 0.2249
    lymp.coef <- -1.3255
    plt.coef <- -0.0052
    crea.coef <- 0.0015
    intercept <- -3.5172
    
    ## predictions
    #predicted.values.death <- pull(data[,match('Sannsynlighetdød', names(data))])
    lin.exp <- intercept + age.coef*age + gender.coef*gender + crp.coef*crp 
    + neut.coef*neut + lymp.coef*lymp + plt.coef*plt + crea.coef*crea
    predicted.values.death <- as.vector(1/(1+exp(-lin.exp)))
    
    ## recalibration (method 1)
    obs.out.D <- mean(outcomeDeath)
    mean.pred.D <- mean(predicted.values.death)
    correction.factor.D <- log((obs.out.D/(1-obs.out.D))/(mean.pred.D/(1-mean.pred.D)))
    intercept <- intercept + correction.factor.D
    lin.exp <- intercept + age.coef*age + gender.coef*gender + crp.coef*crp 
    + neut.coef*neut + lymp.coef*lymp + plt.coef*plt + crea.coef*crea
    predicted.values.death <- as.vector(1/(1+exp(-lin.exp)))
    predval.all.death <- cbind(predval.all.death, predicted.values.death)
    
    ## read data - poor outcome
    data <- read_excel(paste(imputation_dir, '/Zhang/validation_Zhang2_fill_',impu.type,'_',b,'.xlsx',sep=''))
    
    ## outcome
    outcomePoor <- pull(data[,match('ICU or death(inhospital):', names(data))])

    ## covariates
    age <- pull(data[,match('Alder_Ny', names(data))])
    gender <- pull(data[,match('sex', names(data))])
    crp <- pull(data[,match('CRP', names(data))])
    neut <- pull(data[,match('neutrophils', names(data))])
    lymp <- pull(data[,match('Lymfocytter', names(data))])
    plt <- pull(data[,match('platelets', names(data))])
    crea <- pull(data[,match('kreatinin', names(data))])
      
    ## model coefficients
    # (Supplementary Table 2, pag.10-11, Supplementary Material of Zhang et al.)
    age.coef <- 0.0281
    gender.coef <- -0.4723
    crp.coef <- 0.0082
    neut.coef <- 0.2575
    lymp.coef <- -1.1842
    plt.coef <- -0.0027
    crea.coef <- 0.0099
    intercept <- -4.1865
    
    ## predictions

    lin.exp <- intercept + age.coef*age + gender.coef*gender + crp.coef*crp 
    + neut.coef*neut + lymp.coef*lymp + plt.coef*plt + crea.coef*crea
    predicted.values.poor <- as.vector(1/(1+exp(-lin.exp)))
    
    ## recalibration (method 1)
    obs.out.P <- mean(outcomePoor)
    mean.pred.P <- mean(predicted.values.poor)
    correction.factor.P <- log((obs.out.P/(1-obs.out.P))/(mean.pred.P/(1-mean.pred.P)))
    intercept <- intercept + correction.factor.P
    lin.exp <- intercept + age.coef*age + gender.coef*gender + crp.coef*crp 
    + neut.coef*neut + lymp.coef*lymp + plt.coef*plt + crea.coef*crea
    predicted.values.poor <- as.vector(1/(1+exp(-lin.exp)))
    predval.all.poor <- cbind(predval.all.poor, predicted.values.poor)
    
    
  }## end loop on imputation tables
  
  # save
  save(predval.all.death, predval.all.poor, outcomeDeath, outcomePoor, file=paste('./Zhang/validation_Zhang_',type,'.RData',sep=''))
  
}## end loop on imputation method




####---------------------------------------
#### validation of all models on OSLO data:
#### actual plotting
####---------------------------------------

nbins <- 7 

which.model <- c('Allenbach', 'Xie', 'Zhang')
which.impu.type <- c('singleKNN')

for(mymodel in which.model){
  
  # model-specific labels 
  switch(match(mymodel, which.model),
         {myOutcomeLabel <- 'death-ICU'},
         {myOutcomeLabel <- 'death'},
         {myOutcomeLabel <- c('death','poor outcome')})
         
  for(type in which.impu.type){
    
    # load data
    load(paste('./',mymodel,'/validation_',mymodel,'_',type,'.RData',sep=''))
    
    # code preparation 
    switch(match(type, which.impu.type),
           {B <- 0; impu.type <- 'KNN'},
           {B <- 0; impu.type <- 'rf'},
           {B <- 9; impu.type <- 'br'},
           {B <- 9; impu.type <- 'gp'})
    
    for(OutcomeType in myOutcomeLabel){
      
      if(mymodel=='Zhang'){
        # select proper outcome type
        switch(match(OutcomeType, myOutcomeLabel),
               {outcome <- outcomeDeath
               predval.all <- predval.all.death
               modelname <- 'Zhang1'},
               {outcome <- outcomePoor
               predval.all <- predval.all.poor
               modelname <- 'Zhang2'})
      }else{modelname <- mymodel}
      
      ## create objects necessary to plot the ROC curve
      roc.all <- mes.all <- NULL
      for(b in 0:B){
        myROC <- rocit(predval.all[,b+1], outcome)
        roc.all <- c(roc.all, list(myROC))
        myMeasures <- measureit(predval.all[,b+1], outcome, measure = c('PPV', 'NPV'))
        mes.all <- c(mes.all, list(myMeasures))
      }
      
      ## PLOT 1: probability poor outcome
      indJa <- which(outcome==0) # ja = survived! => NO poor outcome
      indNei <- which(outcome==1)
      validation_plot_prefix = './validation_plots_all_models/validation_'
      png(paste(validation_plot_prefix,mymodel,'_',type,'_',OutcomeType,'.png',sep=''), height = 600, width = 750)
      plot(0:3, 0:3, pch="", main = paste(modelname, ' validation, OUH cohort',sep=''), ylim=c(min(rowMeans(predval.all), na.rm = TRUE),max(rowMeans(predval.all), na.rm = TRUE)),
           xlab = OutcomeType, ylab = 'probability of poor outcome', axes = FALSE)
      points(rep(1, length(indJa)) + rnorm(length(indJa),sd=.05), rowMeans(predval.all)[indJa])
      points(rep(2, length(indNei))+ rnorm(length(indNei),sd=.05), rowMeans(predval.all)[indNei])
      axis(2)
      axis(1, at=1:2, label=c('no', 'yes'))
      segments(0.8, mean(rowMeans(predval.all)[indJa], na.rm = TRUE),1.2, mean(rowMeans(predval.all)[indJa], na.rm = TRUE), col=2, lwd=2)
      segments(1.8, mean(rowMeans(predval.all)[indNei], na.rm = TRUE),2.2, mean(rowMeans(predval.all)[indNei], na.rm = TRUE), col=2, lwd=2)
      dev.off()
      
      if(B==0){
        
        ## PLOT 2: validation plot
        # based on fixed bins
        prob.bins.fixed <- seq(0,1,len=nbins+1)
        # based on quantiles
        prob.bins.quantile <- quantile(predval.all[,1], seq(0,1,len=nbins+1))
        for(whichplot in c('fixed', 'quantile')){# two types of calibration plot
          switch(match(whichplot, c('fixed', 'quantile')), {prob.bins <- prob.bins.fixed}, {prob.bins <- prob.bins.quantile})
          xplot <- yplot <- dimEachBin <- NULL
          for(i in 1:nbins){
            Lbound <- prob.bins[i]
            Ubound <- prob.bins[i+1]
            indBind <- which((predval.all[,1]>=Lbound) & (predval.all[,1]<Ubound))
            if(i==nbins)indBind <- c(indBind, which(predval.all[,1]==1))
            dimEachBin <- c(dimEachBin, length(indBind))
            xplot <- c(xplot, mean(predval.all[indBind,1]))
            yplot <- c(yplot, mean(outcome[indBind]))
          }
          mycalibration <- lm(yplot ~ xplot)
          png(paste(validation_plot_prefix,mymodel,'_',type,'_',OutcomeType,'_calibration_',whichplot,'.png',sep=''), height = 400, width = 500)
          plot(xplot, yplot, pch=16, type='b', xlim = c(0,1), ylim = c(-0.1,1.1), xlab = '', 
               ylab = '', cex.axis=2.0)
          abline(a=0, b=1, col="gray50")
          abline(a = mycalibration$coef[1], b = mycalibration$coef[2], lty=2)
          pxplot = seq(0,1,by = 0.05)
          confInt = predict(mycalibration, newdata=data.frame(xplot=pxplot), interval="confidence", level = 0.95)
          lines(pxplot, confInt[,2], col="blue", lty=2)
          lines(pxplot, confInt[,3], col="blue", lty=2)
          if(whichplot=='fixed')text(xplot, yplot+.07, label = as.character(dimEachBin), col=2)
          dev.off()
        }# end loop on the calibration plot type
        
        ## PLOT 3: ROC curve plot
        png(paste(validation_plot_prefix,mymodel,'_',type,'_',OutcomeType,'_ROC.png',sep=''), height = 550, width = 650)
        plot(roc.all[[1]], legend = FALSE, values = TRUE)
        dev.off()
        
        # ROC table
        plotval <- plot(roc.all[[1]], legend = FALSE, values = TRUE)
        dev.off()
        # elements in Alvaro's table:
        # AUC, n_samp, n_feat, tp, tn, fp, fn, sens, spec, ppv, npv, calibration (slope/intercept) 
        youden.cut <- as.numeric(plotval$`optimal Youden Index point`[4])
        youden.ind <- which(roc.all[[1]]$Cutoff == youden.cut)
        youden.ind.m <- which(mes.all[[1]]$Cutoff == youden.cut)
        tab <- c(plotval$AUC,
                 ciAUC(roc.all[[1]])$low,
                 ciAUC(roc.all[[1]])$upp,
                 roc.all[[1]]$pos_count+roc.all[[1]]$neg_count,
                 4,
                 mes.all[[1]]$TP[youden.ind.m],
                 mes.all[[1]]$TN[youden.ind.m],
                 mes.all[[1]]$FP[youden.ind.m],
                 mes.all[[1]]$FN[youden.ind.m],
                 roc.all[[1]]$TPR[youden.ind],
                 1-roc.all[[1]]$FPR[youden.ind],
                 mes.all[[1]]$PPV[youden.ind.m],
                 mes.all[[1]]$NPV[youden.ind.m],
                 mycalibration$coef[1],
                 confint(mycalibration)[1,],
                 mycalibration$coef[2],
                 confint(mycalibration)[2,])
        tab.all <- rbind(tab)
        tab.all <- as.data.frame(tab.all)
        names(tab.all) <- c('AUC', 'CI_AUC_low', 'CI_AUC_upp', 'n_samp', 'n_feat', 'tp', 'tn', 'fp', 'fn', 
                            'sens', 'spec', 'ppv', 'npv', 'int', 'CI_int_low', 'CI_int_upp', 
                            'slope', 'CI_slope_low',  'CI_slope_upp')
        
        
      }else{
        
        ## PLOT 2: validation plot (WHEN MULTIPLE IMPUTATION)
        # based on fixed bins
        prob.bins.fixed <- seq(0,1,len=nbins+1)
        for(whichplot in c('fixed', 'quantile')){# two types of calibration plot
          Xplot <- Yplot <- XYcoeff <- XYcoeffCIint <- XYcoeffCIslope <- NULL
          for(b in 0:B){
            # based on quantiles
            prob.bins.quantile <- quantile(predval.all[,b+1], seq(0,1,len=nbins+1))
            switch(match(whichplot, c('fixed', 'quantile')), {prob.bins <- prob.bins.fixed}, {prob.bins <- prob.bins.quantile})
            # compute single-dataset calibration
            xplot <- yplot <- dimEachBin <- NULL
            for(i in 1:nbins){
              Lbound <- prob.bins[i]
              Ubound <- prob.bins[i+1]
              indBind <- which((predval.all[,b+1]>=Lbound) & (predval.all[,b+1]<Ubound))
              if(i==nbins)indBind <- c(indBind, which(predval.all[,b+1]==1))
              dimEachBin <- c(dimEachBin, length(indBind))
              xplot <- c(xplot, mean(predval.all[indBind,b+1]))
              yplot <- c(yplot, mean(outcome[indBind]))
            }
            mycalibration <- lm(yplot ~ xplot)
            Xplot <- cbind(Xplot, xplot)
            Yplot <- cbind(Yplot, yplot)
            XYcoeff <- cbind(XYcoeff, mycalibration$coeff)
            XYcoeffCIint <- rbind(XYcoeffCIint, confint(mycalibration)[1,] )
            XYcoeffCIslope <- rbind(XYcoeffCIslope, confint(mycalibration)[2,] )
          }
          # calibration plot (all together)
          mycol <- rainbow(B+3)[-c(1,B+3)]
          png(paste(validation_plot_prefix,mymodel,'_',type,'_',OutcomeType,'_calibration_',whichplot,'.png',sep=''), height = 400, width = 500)
          plot(xplot, yplot, pch="", xlim = c(0,1), ylim = c(-0.1,1.1), xlab = 'predicted probabilities', 
               ylab = 'expected true probabilities', main = paste(modelname,' calibration plot, OUH cohort',sep=''))
          for(b in 0:B) points(Xplot[,b+1], Yplot[,b+1], pch = 16, col = mycol[b+1], type = 'b')
          abline(a=0, b=1, col = 'gray50')
          for(b in 0:B) abline(a = XYcoeff[1, b+1], b = XYcoeff[2, b+1], lty=2, col = mycol[b+1])
          legend("topleft", col = mycol, paste('Multiple imputation ROC ',0:B,sep=''), lwd = 2, bty = 'n', cex = .6)
          legendtext <- c(paste('Mean intercept = ', round(rowMeans(XYcoeff)[1],2),sep=''),
                          paste('Mean slope = ', round(rowMeans(XYcoeff)[2],2),sep=''))
          legend("bottomright", col = c(1,1), legendtext, lty = 2, bty = 'n')
          dev.off()
        }# end loop on the calibration plot type
        
        ## PLOT 3: ROC curve plot (WHEN MULTIPLE IMPUTATION)
        png(paste(validation_plot_prefix,mymodel,'_',type,'_',OutcomeType,'_ROC.png',sep=''), height = 550, width = 650)
        plot(roc.all[[1]], legend = FALSE, values = FALSE, col=c(mycol[1], "gray50"), YIndex = FALSE)
        for(b in 1:B)lines(roc.all[[b+1]]$TPR~roc.all[[b+1]]$FPR, col = mycol[b+1])
        legend("bottomright", col = mycol, paste('Multiple imputation ROC ',0:B,sep=''), lwd = 2, bty = 'n', cex = .6)
        dev.off()
        
        # ROC table
        # AUC, n_samp, n_feat, tp, tn, fp, fn, sens, spec, ppv, npv
        tab.all <- NULL
        for(b in 0:B){
          
          plotval <- plot(roc.all[[b+1]], legend = FALSE, values = TRUE)
          dev.off()
          youden.cut <- as.numeric(plotval$`optimal Youden Index point`[4])
          youden.ind <- which(roc.all[[b+1]]$Cutoff == youden.cut)
          youden.ind.m <- which(mes.all[[b+1]]$Cutoff == youden.cut)
          tab <- c(plotval$AUC,
                   ciAUC(roc.all[[b+1]])$low,
                   ciAUC(roc.all[[b+1]])$upp,
                   roc.all[[b+1]]$pos_count+roc.all[[b+1]]$neg_count,
                   4,
                   mes.all[[b+1]]$TP[youden.ind.m],
                   mes.all[[b+1]]$TN[youden.ind.m],
                   mes.all[[b+1]]$FP[youden.ind.m],
                   mes.all[[b+1]]$FN[youden.ind.m],
                   roc.all[[b+1]]$TPR[youden.ind],
                   1-roc.all[[b+1]]$FPR[youden.ind],
                   mes.all[[b+1]]$PPV[youden.ind.m],
                   mes.all[[b+1]]$NPV[youden.ind.m],
                   XYcoeff[1, b+1],
                   XYcoeffCIint[b+1,],
                   XYcoeff[2, b+1],
                   XYcoeffCIslope[b+1,])
          tab.all <- rbind(tab.all, tab)
          
        }
        tab.all <- as.data.frame(tab.all)
        names(tab.all) <- c('AUC', 'CI_AUC_low', 'CI_AUC_upp', 'n_samp', 'n_feat', 'tp', 'tn', 'fp', 'fn', 
                            'sens', 'spec', 'ppv', 'npv', 'int', 'CI_int_low',  'CI_int_upp',
                            'slope',  'CI_slope_low', 'CI_slope_upp')
        rownames(tab.all) <- 0:B
        
      }# end if (yes/no multiple imputation)
      write.table(tab.all, file = paste(validation_plot_prefix,mymodel,'_',type,'_',OutcomeType,'_table.txt',sep=''))
            
    }# end loop on the different outcomes
  }# end loop on the imputation type
}# end loop on the model
