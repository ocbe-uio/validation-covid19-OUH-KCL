### validation_postprocessing_plots.R
###
### Author: Valeria Vitelli
### Contributers: Andrew Reiner, Alvaro Köhn-Luque
###
### Plots results of validation on King's College Hospital cohort
###
###
### Working directory should contain validation_results.Rdata, the output from
### validation_code.R
###
### This script writes output to  subdirectories 'results' and 'results/plots',
### and expects them to exist.


### rm(list=objects())

library(readxl)
library(dplyr)
library(ROCit)


### load data
load('./validation_results.Rdata')

### preliminaries
nbins <- 20 
# alternative models
which.model <- c('Allenbach', 'Xie', 'Zhang')
# alternative imputation methods
which.impu.type <- c('singleKNN')

### plotting and postprocessing loop
for(mymodel in which.model){
  
  # model-specific labels 
  switch(match(mymodel, which.model),
         {myOutcomeLabel <- 'death-ICU'},
         {myOutcomeLabel <- 'death'},
         {myOutcomeLabel <- c('death','poor outcome')})
  
  # load the correct data for the model
  results <- get(paste('res',mymodel,sep=''))
  
  for(type in which.impu.type){
    
    # code preparation for the different imputations (add options if needed!!)
    switch(match(type, which.impu.type),
           {B <- 0; impu.type <- 'KNN'})
    
    # loop on the outcome (Zhang has 2 outcomes)
    for(OutcomeType in myOutcomeLabel){
      
      if(mymodel=='Zhang'){
        # select proper outcome type
        switch(match(OutcomeType, myOutcomeLabel),
               {outcome <- results$out[,1]
               predval.all <- results$pred[,1, drop=FALSE]
               modelname <- 'Zhang1'},
               {outcome <- results$out[,2]
               predval.all <- results$pred[,2, drop=FALSE]
               modelname <- 'Zhang2'})
      }else{
        outcome <- results$out
        predval.all <- cbind(results$pred)
        modelname <- mymodel
      }
      
      # create objects necessary to plot the ROC curve
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
      png(paste('./results/plots/validation_',mymodel,'_',type,'_',OutcomeType,'.png',sep=''), height = 600, width = 750)
      plot(0:3, 0:3, pch="", main = paste(modelname, ' validation, KCH cohort',sep=''), ylim=c(min(rowMeans(predval.all), na.rm = TRUE),max(rowMeans(predval.all), na.rm = TRUE)),
           xlab = OutcomeType, ylab = paste('probability of ',OutcomeType,sep=''), axes = FALSE)
      points(rep(1, length(indJa)) + rnorm(length(indJa),sd=.05), rowMeans(predval.all)[indJa])
      points(rep(2, length(indNei))+ rnorm(length(indNei),sd=.05), rowMeans(predval.all)[indNei])
      axis(2)
      axis(1, at=1:2, label=c('no', 'yes'))
      segments(0.8, mean(rowMeans(predval.all)[indJa], na.rm = TRUE),1.2, mean(rowMeans(predval.all)[indJa], na.rm = TRUE), col=2, lwd=2)
      segments(1.8, mean(rowMeans(predval.all)[indNei], na.rm = TRUE),2.2, mean(rowMeans(predval.all)[indNei], na.rm = TRUE), col=2, lwd=2)
      dev.off()
      
      if(B==0){# currently the only option (only KNN single imputation).. add "else" if needed!!
        
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
          png(paste('./results/plots/validation_',mymodel,'_',type,'_',OutcomeType,'_calibration_',whichplot,'_KCH.png',sep=''), height = 400, width = 500)
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
        png(paste('./results/plots/validation_',mymodel,'_',type,'_',OutcomeType,'_ROC.png',sep=''), height = 550, width = 650)
        plot(roc.all[[1]], legend = FALSE, values = TRUE)
        dev.off()
        
        # ROC table
        plotval <- plot(roc.all[[1]], legend = FALSE, values = TRUE)
        dev.off()
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
        
      }
      write.table(tab.all, file = paste('./results/validation_',mymodel,'_',type,'_',OutcomeType,'_table.txt',sep=''))
      
    }# end loop on the different outcomes
  }# end loop on the different imputation methods
}# end loop on the different models

