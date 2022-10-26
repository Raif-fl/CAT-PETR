#   Copyright (C) 2004-2005 Suman Sundaresh, Computer Science Department
#   University of California, Irvine (suman@uci.edu)
#   Institute for Genomics and Bioinformatics
#   Parts of this code obtained or adapted from "hdarray" (Tony D. Long)

#   This code is free for academic, non-commercial, research use only. If you use this library, please
#   cite the reference (Baldi and Long, 2001) - details are on the Cyber-T
#   Help webpage. Also, please keep all headers in this code as is.
#   For commercial licenses, please contact Pierre Baldi (pfbaldi@uci.edu).


#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#   bayesreg: v2.0beta

ALMOST_ZERO = 10E-14;
################################################################################
##   bayesT
##   
##       Adata - Actual data frame (ncol = numC + numE) contiguous controls then experimentals       
##       ppde = FDR test.., fits a beta and uniform mixture dist to the resulting p-values.
##           BETAFIT - param for ppde
##       bayes = do bayesian est of variance
##       winSize = how big around sorted dpoints for prior calc, has to be odd.
##       conf = degrees of freedom for each population prior calculation.
##              (how many dpoints do we believe make up the prior?)
##       doMulttest - allow BH and Bonferroni adjusted p-vals to be 
##                   calculated right here.
##       bayesIntC - use similar means (instead of similar variances) to calculate bayesian variance 
##                   for the Controls
##       bayesIntE - use similar means (instead of similar variances) to calculate bayesian variance
##                   for the Experimentals
################################################################################
bayesT <- function (aData, numC, numE, ppde=TRUE, betaFit=1, bayes=TRUE, winSize=101, conf=10, 
                    doMulttest=FALSE, bayesIntC=FALSE, bayesIntE=FALSE){
  if ((ceiling((winSize-1)/2))!=((winSize-1)/2))
    stop("ERROR: winSize must be an odd number.")
  
  numGene<- nrow(aData)
  
  ## compute number of valid entries for each gene
  nC <- apply(as.matrix(aData[, 1:numC]), 1, function(x) sum(!is.na(x)))
  nE<- apply(as.matrix(aData[, (numC+1):(numC+numE)]), 1, function(x) sum(!is.na(x)))
  
  ## compute means for valid entries
  meanC<- apply(as.matrix(aData[, 1:numC]), 1, function(x) 
                if (sum(!is.na(x)))
                mean(x[!is.na(x)])
                else NA)
  meanE<- apply(as.matrix(aData[, (numC+1):(numC+numE)]), 1, function(x) 
                if (sum(!is.na(x)))
                mean(x[!is.na(x)])
                else NA)
  stdC<- apply(as.matrix(aData[, 1:numC]), 1, function(x) 
               if (sum(!is.na(x)) > 1)
               sqrt(var(x[!is.na(x)]))
               else NA)
  stdE<- apply(as.matrix(aData[, (numC+1):(numC+numE)]), 1, function(x) 
               if (sum(!is.na(x)) > 1)
               sqrt(var(x[!is.na(x)]))
               else NA)

  ## Calc the statistics in a Bayesian or Classical way
  if (bayes){
    ## Do something a little different if we only have a single replicate
    if (numC == 1) {bayesIntC <- TRUE};
    if (numE == 1) {bayesIntE <- TRUE};
    
    ## Here we calculate the average sds over windows of probes.
    rasdC <- rep(NA, numGene)
    if (!bayesIntC){
      ## basic case, order by means, averge the sds in a windowq.
      temp <- runavg(stdC[!is.na(stdC)][order(meanC[!is.na(stdC)])],((winSize-1)/2))	
      temp <- temp[rank(meanC[!is.na(stdC)])]
      rasdC[!is.na(stdC)]<- temp
    } else {
      ## use the actual intensities ranked by the means
      #cat('About to calc the bayes Variance for control using intensities.\n')
      intMat <- as.matrix(aData[!is.na(meanC), 1:numC]);
      intMat <- intMat[order(meanC[!is.na(meanC)]), ];
      temp <- runavgPool(intMat, winSize)
      temp <- temp[rank(meanC[!is.na(meanC)])]
      rasdC[!is.na(meanC)]<- temp
    }
    
    ## Same here, check for a single replicate in the experiment.
    rasdE <- rep(NA, numGene)
    if (!bayesIntE){
      temp <- runavg(stdE[!is.na(stdE)][order(meanE[!is.na(stdE)])],((winSize-1)/2))	
      temp <- temp[rank(meanE[!is.na(stdE)])]	
      rasdE[!is.na(stdE)]<- temp
    } else {
      #cat('About to calc the bayes Variance for experimental using intensities\n')
      intMat<- as.matrix(aData[!is.na(meanE), (numC+1):(numC+numE)]);
      intMat <- intMat[order(meanE[!is.na(meanE)]), ];
      temp <- runavgPool(intMat, winSize)	
      temp <- temp[rank(meanE[!is.na(meanE)])]
      rasdE[!is.na(meanE)]<- temp
    }
    
    ## Before computing the Bayes SD, go through and set StdC to almost (zero)  whereever is.na and rasdC is not NA
    ##   Vice versa for the expt data
    forBayesStdC <- stdC
    forBayesStdE <- stdE
    
    forBayesStdC[is.na(stdC) & !is.na(rasdC)] <- ALMOST_ZERO;
    forBayesStdE[is.na(stdE) & !is.na(rasdE)] <- ALMOST_ZERO;
    
    ## Set the ALMOST_ZERO for all rasdC and rasdE at zero as well.
    rasdC[rasdC == 0] <- ALMOST_ZERO;
    rasdE[rasdE == 0] <- ALMOST_ZERO; 
    
    ## compute bayes sd
    if (all(conf + nC == 2) || all(conf + nE == 2)){
      stop('ERROR: Zero degrees of freedom, increase the confidence parameter.');
    }
    bayesSDC<- sqrt((conf * rasdC^2 + (nC - 1) * forBayesStdC^2)/(conf + nC - 2))
    bayesSDE<- sqrt((conf * rasdE^2 + (nE - 1) * forBayesStdE^2)/(conf + nE - 2))
    		
    ##For here, they were using the regular Dfs to calculate the original T value,
    ##   But if we have only a single replicate, then we need to tweak this.., 
    ##   Note with a single rep, we are in effect increasing conf by 1.
    newNumC <- nC
    if (numC < 2){
      newNumC <- 2;
    } 
    newNumE <- nE;
    if (numE < 2){
      newNumE <- 2;
    }
    
    sumStats <- cbind(newNumC,newNumE,meanC,meanE,bayesSDC,bayesSDE)
    ttest <- t(apply(sumStats, 1, function(x) tstat(x)))
    colnames(ttest)=c("bayesT","bayesDF","varRatio")
    ##change Bayes degree of freedom to reflect pseudo-counts
    ttest[,2] <- ttest[,2] + 2 * conf - 2		    
    
    ##Subtract out the fake df for each single rep case:
    toSub <- sum(c(numC ==1, numE == 1 ));
    ttest[, 2] <- ttest[, 2] - toSub;
  } else {
    ## A regular t-test
    sumStats <- cbind(nC, nE, meanC, meanE, stdC, stdE)
    ttest <- t(apply(sumStats, 1, function(x) tstat(x)))
    colnames(ttest) = c("T", "DF", "varRatio")
  }
  
  ## Pvals amd fold changes
  pVal <- 1 - pf(ttest[,1]^2, 1, ttest[,2])
  fold <- - (meanC/meanE) * ((meanE/meanC) < 1) + (meanE/meanC) * ((meanC/meanE) < 1)
  
  ## Then the final data.frame 
  if (bayes){
    objBayes<- cbind(aData, nC, nE, meanC, meanE, stdC, stdE, fold, rasdC, rasdE, bayesSDC, bayesSDE,ttest, 
                     pVal)
  }
  else{
    objBayes<- cbind(aData, nC, nE, meanC, meanE, stdC, stdE, fold, ttest, pVal)
  }

  ##posterior probability of differential expression (ppde) test if wanted
  if (ppde){
    ppdeVal <- ppdeMix(as.matrix(pVal), betaFit)
    colnames(ppdeVal) = c("cum.ppde.p", "ppde.p", "ROC.x", "ROC.y")
    objBayes <- cbind(objBayes, ppdeVal)
  }

  ## Mult testing correction
  if (doMulttest){
    ##Bind them, last two cols of adjp list in original order
    objBayes <- data.frame(objBayes, runMulttest(pVal))
  }
  
  return(objBayes)
}

################################################################################
## Simple wrapper to load and then run vsn
################################################################################
runVsn <- function(data, ...){
  .libPaths("/home/baldig/shared_libraries/centos64/pkgs/R/2.15.1_fix/lib64/R/library")
  suppressMessages(require(vsn))
  #cat('data in VSN: \n', colnames(data), '\n', as.vector(unlist(data[1,])), '\n',
  #    as.vector(unlist(data[2,])), '\n', as.vector(unlist(data[3, ])), '\n');
  tryCatch(data.frame(justvsn(as.matrix(data), verbose=FALSE, minDataPointsPerStratum=10, ...)),
           error=function(e){
             libs_a<-.libPaths()
             writeError('Error in running vsn.  Possible collinearity.  Check your data.')
             stop(e)
           }
         )
}

################################################################################
## Simple function to do Bonferroni and BH multiple testing correction.
################################################################################
runMulttest <- function(pvals){
  .libPaths("/home/baldig/shared_libraries/centos64/pkgs/R/2.15.1_fix/lib64/R/library")
  suppressMessages(require(multtest));
  adjPObj <- mt.rawp2adjp(pvals, proc=c("Bonferroni", "BH"))
  ##Bind them, last two cols of adjp list in original order
  adjPObj$adjp[order(adjPObj$index), -1]
}

################################################################################
## bayesT.pair
## 
## NOTE: aData needs to be log-transformed
## estExp is the last column of aData and it represents
## estimated expression for that gene 
## A good value for 'estExp; would be the mean of the log of the
## "extimated expression level" over both treatments (i.e., control and
## experimentals) and replicates, where the estimated expression level for each
## gene/treatment/replicate is given as a fraction of total expression over
## all genes for that treatment/replicate.
## 
## doMulttest - do Bonf and BH correction?
## bayesIntR - pool ratios for the variance, instead of pooling variances
## 
################################################################################
bayesT.pair <-  function (aData, numR, ppde=TRUE, betaFit=1, bayes=TRUE, winSize=101, conf=10,
                          doMulttest=FALSE, bayesIntR=FALSE){
  if ((ceiling((winSize-1)/2))!=((winSize-1)/2))
    stop("ERROR: winSize must be an odd number.")
  
  estExp <- numR + 1
  numGene <- nrow(aData)
  
  ##compute number of valid entries for each gene
  nR <- apply(as.matrix(aData[, 1:numR]), 1, function(x) sum(!is.na(x)))

  ##compute means for valid entries
  meanR <- apply(as.matrix(aData[, 1:numR]),
                1, function(x) {
                  if (sum(!is.na(x)))
                    mean(x[!is.na(x)])
                  else NA}
                )
  stdR <- apply(as.matrix(aData[, 1:numR]),
               1, function(x){ 
                 if (sum(!is.na(x)) > 1)
                   sqrt(var(x[!is.na(x)]))
                 else NA}
               )

  ## Calc the statistic in a bayesian or classical manner
  if (bayes){
    rasdR <- rep(NA, numGene)
    index.col <- aData[,estExp]
    
    ## If we are using intensities, do things a little differently
    if (numR == 1){bayesIntR <- TRUE}
    
    if (!bayesIntR){
      ## Calculate averages over sd windows
      temp <- runavg(stdR[!is.na(stdR)][order(index.col[!is.na(stdR)])],((winSize-1)/2))	
      temp <- temp[rank(index.col[!is.na(stdR)])]	
      rasdR[!is.na(stdR)]<- temp
    } else {
      ## Calculate averages of sds over windows of intensities
      #cat('About to calc bayes Variance for ratios using intensities.')
      intMat <- as.matrix(aData[!is.na(meanR), 1:numR])
      intMat <- intMat[order(meanR[!is.na(meanR)]), ]
      temp <- runavgPool(intMat, winSize)
      temp <- temp[rank(meanR[!is.na(meanR)])]
      rasdR[!is.na(meanR)] <- temp
    }

    forBayesStdR <- stdR
    forBayesStdR[is.na(stdR)] <- ALMOST_ZERO
    
    ##compute bayes sd
    bayesSD<- sqrt((conf * rasdR^2 + (nR - 1) * forBayesStdR^2)/(conf + nR - 2))
    ttest<- sqrt(nR) * (meanR/bayesSD)
    
    ##change Bayes degree of freedom to reflect pseudo-counts
    bayesDF <- nR + conf - 2
    pVal <- 1 - pf(ttest^2, 1, bayesDF)
  }
  else{
    ## Calc the basic paired tstats
    ttest <- sqrt(nR) * (meanR/stdR)
    DF <- nR - 1
    pVal <- 1 - pf(ttest^2, 1, DF)
  }

  ## Put together the final version of the output files 
  if (bayes){
    objBayes<- data.frame(aData, nR, meanR, stdR, rasdR, bayesSD, ttest, bayesDF,pVal)
  }
  else{
    objBayes<- data.frame(aData, nR, meanR, stdR, ttest, DF, pVal)
  }
  
  ##posterior probability of differential expression (ppde) test if wanted
  if(ppde){
    ppdeVal <- ppdeMix(as.matrix(pVal),betaFit)
    colnames(ppdeVal)=c("cum.ppde.p","ppde.p","ROC.x","ROC.y")
    objBayes<-data.frame(objBayes, ppdeVal)
  }

  ## Multtest if wanted
  if (doMulttest){
    objBayes <- data.frame(objBayes, runMulttest(pVal))
  }
  
  return(objBayes)
}

################################################################################
## Helper functions for the tstat calculations
################################################################################
tstat <- function (sumStats)
{
  ##number of replicates, means and std deviations
  n1=sumStats[1];  n2=sumStats[2]
  m1=sumStats[3];  m2=sumStats[4]
  sd1=sumStats[5]; sd2=sumStats[6]
  ##  do the two sample t-test
  tt <- -(m1 - m2)/sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 +
                                                                n2 - 2)) * ((n1 + n2)/(n1 * n2)))
  dft <- n1 + n2 - 2
  rvar <- max((sd1^2)/(sd2^2), (sd2^2)/(sd1^2))
  as.vector(c(tt, dft, rvar))
}


runavg <- function(x,k=1) {
  ## x is the input vector
  ## n is the length of the vector
  ## k is the size of the running average window, the window 
  ##   includes the data point at that position plus k points 
  ##   on either side.  For example k = 3 would be a sliding
  ##   window of size {3+1+3}= 7 centered on each point
  ## r is the output vector
  if (k <= 0)
    stop("Error: Window size (k) must be greater than zero")
  x<-as.array(x);
  n <- length(x)
  if (k >= n){
    msg <- 'Window Size must be less than number of data points!'
    writeError(msg)
    stop(msg)
  }
  
  r<-0
  t<-0
  
  for (j in 0:k) {
    t<-t+x[j+1]
  }
  
  l=0
  u=k
  
  for (i in 0:(n-1)) {
    ## the use of u-l-1 instead of k handles end effects 
    r[i+1]= (t/(u-l+1))
    ## update the current running total for the next position 
    if (i > k-1) {
      t=t-x[l+1]
      l=l+1
    }
    if (u < n-1) {
      u=u+1			
      t=t+x[u+1]
    }
  }
  
  return(as.array(r))
}

runavgPool <- function(x,k=1) {
  ## x is the input vector
  ## n is the length of the vector
  ## k is the size of the running average window, the window 
  ##   includes the data point at that position plus k points 
  ##   on either side.  For example k = 3 would be a sliding
  ##   window of size {3+1+3}= 7 centered on each point
  ## r is the output vector
  if (k <= 0)
    stop("Error: Window size (k) must be greater than zero")
  
  x <- as.matrix(x);
  n <- dim(x)[1]
  
  r<-0
  t<-0
  
  side <- (k-1)/2
  
  for (j in 1:side) {
    flatArr <- NULL
    for(i in 1:(j+side))
      flatArr <- c(flatArr,x[i,])
    
    flatArr<-as.array(flatArr)
    
    r[j]<-sd(flatArr)		
  }
  for (j in (side+1):(n-(side+1))) {
    flatArr <- NULL
    for (i in (j-side):(j+side)) 
      flatArr <- c(flatArr,x[i,])
    flatArr<-as.array(flatArr)	
    r[j]<-sd(flatArr)		
    
  }
  for (j in (n-side):n) {
    flatArr <- NULL
    for(i in (j-side):n)
      flatArr <- c(flatArr,x[i,])
    flatArr<-as.array(flatArr)	
    r[j]<-sd(flatArr)		
  }
  return(as.array(r))
}

secMin<-function(a){
  ##returns the next minimum after the overall min
  i=0
  a<-sort(a)
  secondMin=a[1]
  while (secondMin==a[1]){
    i=i+1
    secondMin=a[i]
  }
  return(secondMin)
}

secMax<-function(a){
  ##returns the next maximum after the overall max
  i=0
  a<-sort(a,decreasing=TRUE)
  secondMax=a[1]
  if (all(a == secondMax)){
    return(NA)
  }
  while (secondMax==a[1]){
    i=i+1
    secondMax=a[i]
  }
  return(secondMax)
}


adjPVal<-function(a){
  ##makes sure that p-values are not zero but some very small value (factor of 10 smaller than second min)
  ##also if p-values are 1, then they are made a very small number less than 1 but more than all others
  if(min(a)==0){
    if (all(a == 0)){a[a==0] = 0.0000000001;  return(a)}
    adj<-secMin(a)/10
    a[a==0]=adj
  }
  if(max(a)==1){
    if (all(a == 1)){a[a==1] = 0.999999999;  return(a)}
    
    adj<-secMax(a) + ((1-secMax(a))*9/10)
    a[a==1]=adj
  }
  return(a)
}


################################################################################
## cyberTPlots
## Make some simple plots for the webserver
## resData - has normalized data, pVals, ROC vals
## rawData - has original data sent in.
## Make the following plots:
## Mean Control Vs Mean Experimental, raw and norm
## Mean vs Variance, raw and norm, mean and control
## p-value associated with rawScores (fold change (or mean diff))
## ROC Plot (if given)
## Some plots only make sense for CE or paired, so take what we can.
## Assumes that the dataframes are lined up! And dfs are in format:
##  Controls, Experiments, Results
## And resData - has normalized data
##
## 4-8-12 - Add in smoothing and outlier removal options
## 1-25-12 - Make so write out each in different png files for display on web
##
################################################################################
cyberTPlots <- function(resData, rawData, numC=0, numE=0, numR=0, repPerCond=NULL, dataType='CE', normType='None',
                        doSmoothing=TRUE, doOutlierRemoval=TRUE, outputFmt='png', doTitle=TRUE){
  .libPaths("/home/baldig/shared_libraries/centos64/pkgs/R/2.15.1_fix/lib64/R/library")
  library(lattice)
  library(geneplotter)

  #fileNameFmt <- sprintf('%%s.%s', outputFmt)
  fileNameFmt <- outputFmt
  graphicsFunc <- png
  if (outputFmt == 'ps'){
    graphicsFunc <- postscript
  }

  mainFunc <- function(title){
    if (doTitle){
      return(title)
    }
    NULL
  }
  
  ## Setup functions for smoothing and outlier removal
  smoothingPanelFunc <- function(...){panel.smoothScatter(nbin=128, ...)}
  outlierRows <- function(x){
    qnt <- quantile(x, probs=c(.25, .75), na.rm=TRUE)
    h <- 2.0 * IQR(x, na.rm=TRUE)
    is.na(x) | is.nan(x) | (x < (qnt[1] - h)) | (x > (qnt[2] + h)) 
  }
  outlierRows.2d <- function(x,y){
    outlierRows(x) | outlierRows(y)
  }
  outlierRows.splom <- function(d){
    alld <- as.vector(unlist(d))
    qnt <- quantile(alld, probs=c(.25, .75), na.rm=TRUE)
    h <- 2.0 * IQR(alld, na.rm=TRUE)
    outDF <- is.na(d) | (d < qnt[1] - h) | (d > qnt[2] + h) 
    apply(outDF, 1, any)
  }

  panelFunc <- panel.xyplot #lattice.getOption("panel.xyplot")
  splomPanelFunc <- panel.splom #lattice.getOption("panel.splom")
  outlierFunc <- function(...){FALSE}
  outlierSplomFunc <- function(...){FALSE}
  if (doSmoothing){
    panelFunc <- smoothingPanelFunc
    splomPanelFunc <- smoothingPanelFunc
  }
  if (doOutlierRemoval){
    outlierFunc <- outlierRows.2d
    outlierSplomFunc <- outlierRows.splom
  }

  getPanelFunc <- function(){
    return(panelFunc)
  }
  
  cCols <- c(1:numC)
  eCols <- c((numC+1):(numC+numE))
  rCols <- c(1:numR)
  
  bayes <- FALSE
  if (dataType == 'CE'){
    ## Setup the data to plot
    if (normType == 'None'){
      toPlot <- data.frame(type=factor(rep(c('Cond1', 'Cond2'), each=nrow(rawData))),
                           rawmean=c(resData$meanC, resData$meanE),
                           rawsd=c(resData$stdC, resData$stdE),
                           normmean=c(resData$meanC, resData$meanE),
                           normsd=c(resData$stdC, resData$stdE)
                           )
      rawMeans <- data.frame(C1=resData$meanC, C2=resData$meanE) 
    } else if (normType != 'None'){
      toPlot <- data.frame(type=factor(rep(c('Cond1', 'Cond2'), each=nrow(rawData))),
                           rawmean=c(apply(rawData[,cCols], 1, mean), apply(rawData[,eCols], 1, mean)),
                           rawsd=c(apply(rawData[,cCols], 1, sd), apply(rawData[,eCols], 1, sd)),
                           normmean=c(resData$meanC, resData$meanE),
                           normsd=c(resData$stdC, resData$stdE))
      rawMeans <- data.frame(C1=apply(rawData[,cCols], 1, mean), C2=apply(rawData[,eCols], 1, mean))
    }
    if (length(grep('bayes', colnames(resData))) > 0){
      toPlot <- data.frame(toPlot,
                           bayessd=c(resData$bayesSDC, resData$bayesSDE))
      bayes <- TRUE;
    }

    ## Start with some plots
    ## Means of different conditions
    graphicsFunc(sprintf('ce_raw_plot.%s', fileNameFmt))
    print(xyplot(C2 ~ C1, rawMeans, xlab = "Mean Conditon 1 Raw Intensity", ylab = "Mean Condition 2 Raw Intensity",
                 main = mainFunc("Raw Intensity Comparison"),
                 panel=panelFunc, subset=!outlierFunc(rawMeans$C1, rawMeans$C2)))
    graphics.off()
    
    if (normType != 'None'){
      normMeans <- data.frame(C1=resData$meanC, C2=resData$meanE)
      graphicsFunc(sprintf('ce_norm_plot.%s', fileNameFmt))
      xlab <- sprintf("Mean Condition 1 %s-transformed Intensity", normType)
      ylab <- sprintf("Mean Condition 2 %s-transformed Intensity", normType)
      print(xyplot(C2 ~ C1, normMeans, 
             xlab = xlab, ylab = ylab,
             main = mainFunc(sprintf("Normalized Intensity Comparison", normType)),
                   panel=panelFunc, subset=!outlierFunc(normMeans$C1, normMeans$C2)))
      graphics.off()
    }

    ## Means vs sds
    ## Only if Norm is not None or no bayes
    ## if norm is None and bayes, then this is included in the bayes plot below.
    if (normType != 'None' || ! bayes)
    {
      graphicsFunc(sprintf('sdmean_raw_plot.%s', fileNameFmt))
      print(xyplot(rawsd ~ rawmean | type, toPlot, xlab = "Mean Raw Intensity",
                   ylab = "Empirical Standard Deviation ",
                   main = mainFunc("Raw Mean/Empirical Variance Dependence"),
                   panel=panelFunc, subset=!outlierFunc(toPlot$rawsd, toPlot$rawmean)))
      graphics.off()
    }

    ## And for the next, only plot if no bayes.
    ## If there is bayes, then we will have included in another plot.
    if (normType != 'None' && !bayes)
    {
      graphicsFunc(sprintf('sdmean_norm_plot.%s', fileNameFmt))
      xlab <- sprintf("Mean %s-transformed Intensity", normType)
      print(xyplot(normsd ~ normmean | type, toPlot,
                   xlab = xlab,
                   ylab = "Empirical Standard Deviation",
                   main = mainFunc(sprintf("%s-transformed Mean/Empirical Variance Dependence", normType)),
                   pch="*", panel=panelFunc, subset=!outlierFunc(toPlot$normsd, toPlot$normmean)))
      graphics.off()
    }
    
    if (bayes){
      ## Plot mean vs bayessd  - Actually, this is in the plot lower down
      #graphicsFunc(sprintf('bayessdmean_norm_plot.%s', fileNameFmt))
      #xlab <- "Mean Raw Intensity"
      #main <- "Raw Mean/Bayes-regularized Variance Dependence"
      #if (normType != 'None'){
      #  xlab <- sprintf("Mean %s-transformed Intensity", normType)
      #  main <- sprintf("%s-transformed Mean/Bayes-regularized Variance Dependence", normType)
      #}
      #print(xyplot(bayessd ~ normmean| type, toPlot, ylab='Bayes-regularized Standard Deviation',
      #             xlab=xlab,
       #             main=main
      #             ))
      #graphics.off()

      ## Plot empirical vs bayes sds
      graphicsFunc(sprintf('bayesvempsd_norm_plot.%s', fileNameFmt))
      main <- mainFunc("Raw Data Regularization Effect")
      if (normType != 'None'){
        main <- sprintf("%s-transformed Data Regularization Effect", normType)
      }
      print(xyplot(bayessd ~ normsd |type, toPlot, xlab='Empirical Standard Deviation',
                   ylab='Bayes-regularized Standard Deviation',
                   main=main ,
                   panel = function(...) {
                     panelFunc(...);
                     panel.abline(a = c(0, 1), lty = 2);
                   }, subset=!outlierFunc(toPlot$bayessd, toPlot$normsd)))
      graphics.off()

      ## Plot of mean vs std devs for norm data, both bayes and empirical
      typeVect <- c(levels(toPlot$type)[toPlot$type], levels(toPlot$type)[toPlot$type])
      ## And then add in a plot where we look at the mean vs std for bayes/non bayes
      nToPlot <- data.frame(mean=c(toPlot$normmean, toPlot$normmean),
                            sd=c(toPlot$normsd, toPlot$bayessd),
                            sdtype=rep(c('Empirical', 'Bayes-regularized'), each=nrow(toPlot)),
                            type=typeVect)
      xlab <- "Mean Raw Intensity"
      main <- mainFunc("Effect of regularization - Mean/Variance Dependence")
      if (normType != 'None'){
        xlab <- sprintf("Mean %s-transformed Intensity", normType)
      }
      graphicsFunc(sprintf('bayesandempsd_vs_mean_plot.%s', fileNameFmt))
      print(xyplot(sd~ mean | sdtype + type, nToPlot, xlab=xlab,
                   ylab="Standard Deviation", main=main,
                   panel=panelFunc,
                   subset=!outlierFunc(nToPlot$sd, nToPlot$mean)))
      graphics.off()
      
    }
  }
  if (dataType == 'P'){
    if (normType == 'None'){
      toPlot <- data.frame(rawmean=resData$meanR,
                           rawsd=resData$stdR,
                           normmean=resData$meanR,
                           normsd=resData$stdR
                           )
    } else if (normType != 'None'){
      toPlot <- data.frame(rawmean=apply(rawData[,rCols], 1, mean),
                           rawsd=apply(rawData[,rCols], 1, sd), 
                           normmean=resData$meanR,
                           normsd=resData$stdR)
    }
    if (length(grep('bayes', colnames(resData))) > 0){
      toPlot <- data.frame(toPlot,
                           bayessd=resData$bayesSD)
      bayes <- TRUE;
    }
    if (normType != 'None' || ! bayes)
    {
      graphicsFunc(sprintf('sdmean_raw_plot.%s', fileNameFmt))
      print(xyplot(rawsd ~ rawmean, toPlot, xlab = "Raw Mean Difference", ylab = "Empirical Standard Deviation",
                   main = mainFunc("Raw Mean Difference/Empirical Variance Dependence"), panel=panelFunc,
                   subset=!outlierFunc(toPlot$rawsd, toPlot$rawmean)))
      graphics.off()
    }
    if (normType != 'None' && ! bayes){
      graphicsFunc(sprintf('sdmean_norm_plot.%s', fileNameFmt))
      xlab <- sprintf("Mean %s-transformed Difference", normType)
      print(xyplot(normsd ~ normmean, toPlot, xlab = xlab, ylab = "Empirical Standard Deviation",
                   main = mainFunc(sprintf("%s-transformed Mean Difference/Empirical Variance Dependence",
                     normType)), panel=panelFunc, subset=!outlierFunc(toPlot$normsd, toPlot$normmean)))
      graphics.off()
    }
    if (bayes){
      ## graphicsFunc(sprintf('bayessdmean_norm_plot.%s', fileNameFmt))
      ## xlab <- "Raw Mean Differences"
      ## main <- "Raw Mean Difference/Bayes-regularized Variance Dependence"
      ## if (normType != 'None'){
      ##   xlab <- sprintf("Mean %s-transformed Difference", normType)
      ##   main <- sprintf("%s-transformed Mean/Bayes-regularized Variance Dependence", normType)
      ## }
      ## print(xyplot(bayessd ~ normmean, toPlot, ylab='Bayes-regularized Standard Deviation',
      ##              xlab=xlab,
      ##              main=main
      ##              ))
      ## graphics.off()
      graphicsFunc(sprintf('bayesvempsd_norm_plot.%s', fileNameFmt))
      main <- mainFunc("Raw Data Regularization Effect")
      if (normType != 'None'){
        main <- mainFunc(sprintf("%s-transformed Data Regularization Effect", normType))
      }
      print(xyplot(bayessd ~ normsd , toPlot, xlab='Empirical Standard Deviation',
                   ylab='Bayes-regularized Standard Deviation',
                   main=main,
                   panel = function(...) {
                     panelFunc(...);
                     panel.abline(a = c(0, 1), lty = 2);
                   }, subset=!outlierFunc(toPlot$bayessd, toPlot$normsd)))
      graphics.off()

      ## Combining the above :
      #typeVect <- c(levels(toPlot$type)[toPlot$type], levels(toPlot$type)[toPlot$type])
      ## And then add in a plot where we look at the mean vs std for bayes/non bayes
      nToPlot <- data.frame(mean=c(toPlot$normmean, toPlot$normmean),
                            sd=c(toPlot$normsd, toPlot$bayessd),
                            sdtype=rep(c('Empirical', 'Bayes-regularized'), each=nrow(toPlot))
                            )
                            #type=typeVect)
      xlab <- "Mean Raw Differences"
      main <- mainFunc("Effect of regularization - Mean/Variance Dependence")
      if (normType != 'None'){
        xlab <- sprintf("Mean %s-transformed Differences", normType)
      }
      graphicsFunc(sprintf('bayesandempsd_vs_mean_plot.%s', fileNameFmt))
      print(xyplot(sd~ mean | sdtype , nToPlot, xlab=xlab,
                   ylab="Standard Deviation", main=main,
                   panel=panelFunc, subset=!outlierFunc(nToPlot$sd, nToPlot$mean)))
      graphics.off()
      
    }
  }
  if (dataType == 'ANOVA'){
    ## Setup the data to plot
    extractData <- function(nameFmt){
      res <- c()
      for (i in 1:length(repPerCond)){
        colName <- sprintf(nameFmt, i)
        res <- c(res, resData[[colName]])
      }
      res
    }
    if (normType == 'None'){
      toPlot <- data.frame(type=factor(rep(c(sprintf('C%d', 1:length(repPerCond))), each=nrow(rawData))),
                           rawmean=extractData('Mean%d'),
                           rawsd=extractData('SD%d'))
      toPlot <- data.frame(toPlot,
                           normmean=toPlot$rawmean,
                           normsd=toPlot$rawsd
                           )
      rawMeans <- resData[, grep('Mean', colnames(resData))]
    } else if (normType != 'None'){
      toPlot <- data.frame(type=factor(rep(c(sprintf('C%d', 1:length(repPerCond))), each=nrow(rawData))),
                           normmean=extractData('Mean%d'),
                           normsd=extractData('SD%d'))

      colEnds <- cumsum(repPerCond) 
      colStarts <- colEnds - repPerCond + 1
      ## cat('colStarts: ', colStarts, ', colEnds: ', colEnds, '\n')

      rawMeans <- matrix(NA, nrow=nrow(rawData), ncol=length(repPerCond))
      rawSds <- matrix(NA, nrow=nrow(rawData), ncol=length(repPerCond))

      for (i in 1:length(repPerCond)){
        if (repPerCond[i] == 1){
          rawMeans[, i] <- rawData[,colStarts[i]]
          rawSds[, i] <- 0
        } else {
          rawMeans[, i] <- apply(rawData[,colStarts[i]:colEnds[i]], 1, mean)
          rawSds[, i] <- apply(rawData[,colStarts[i]:colEnds[i]], 1, sd)
        }
      }
      rawMeans <- data.frame(rawMeans)
      colnames(rawMeans) <- sprintf('C%d', 1:length(repPerCond))
      rawSds <- data.frame(rawSds)
      colnames(rawSds) <- sprintf('C%d', 1:length(repPerCond))

      toPlot <- data.frame(toPlot,
                           rawmean=as.vector(unlist(rawMeans)),
                           rawsd=as.vector(unlist(rawSds))
                           )
    }
    if (length(grep('bayes', colnames(resData))) > 0){
      toPlot <- data.frame(toPlot,
                           bayessd=extractData('bayesSD%d'))
      bayes <- TRUE;
    }
    ## Start with some plots
    ## Means of different conditions
    cat('About to print first splom!\n')
    cat('Dim(rawMeans)', dim(rawMeans), '\n')
    toPlotSplom <- rawMeans[!outlierSplomFunc(rawMeans),]
    cat('Dim(toPlot)', dim(toPlot), '\n')
    cat(sprintf('%s', head(toPlot)), '\n')
    graphicsFunc(sprintf('anova_raw_plot.%s', fileNameFmt))
    print(splom(toPlotSplom,
                main = mainFunc("Raw Intensity Comparison"),
                panel=splomPanelFunc))
    graphics.off()
    
    if (normType != 'None'){
      normMeans <- resData[, grep('Mean', colnames(resData))]
      cat('About to print second splom!\n')
      toPlotSplom <- normMeans[!outlierSplomFunc(normMeans),]
      cat(sprintf('%s', head(toPlot)), '\n')
      graphicsFunc(sprintf('anova_norm_plot.%s', fileNameFmt))
      print(splom(toPlotSplom, 
                  main = mainFunc(sprintf("%s-transformed Intensity Comparison", normType)),
                  panel=splomPanelFunc))
      graphics.off()
    }

    ## Means vs sds - Only if Norm or no bayes
    if (normType != 'None' || ! bayes)
    {
      if (!any(is.na(toPlot$rawsd))){
        cat('About to print sdmean_raw_plot.!\n')
        graphicsFunc(sprintf('sdmean_raw_plot.%s', fileNameFmt))
        print(xyplot(rawsd ~ rawmean | type, toPlot, xlab = "Mean Raw Intensity ", ylab = "Empirical Standard Deviation",
                     main = mainFunc("Raw Mean/Empirical Variance Dependence"), panel=panelFunc,
                     subset=!outlierFunc(toPlot$rawsd, toPlot$rawmean)))
        graphics.off()
      }
    }
    if (normType != 'None' && ! bayes){
      if (!any(is.na(toPlot$normsd))){
        cat('About to sdmean_norm_plot.!\n')
        graphicsFunc(sprintf('sdmean_norm_plot.%s', fileNameFmt))
        xlab <- sprintf("Mean %s-transformed Intensity", normType)
        print(xyplot(normsd ~ normmean | type, toPlot,
                     xlab = xlab, ylab = "Empirical Standard Deviation",
                     main = mainFunc(sprintf("%s-transformed Mean/Empirical Variance Dependence",
                       normType)),
                     pch="*", panel=panelFunc,
                     subset=!outlierFunc(toPlot$normsd, toPlot$normmean)))
        graphics.off()
      }
    }
    
    if (bayes){
      ## graphicsFunc(sprintf('bayessdmean_norm_plot.%s', fileNameFmt))
      ## xlab <- "Mean Raw Intensity"
      ## main <- "Raw Mean/Bayes-regularized Dependence"
      ## if (normType != 'None'){
      ##   xlab <- sprintf("Mean %s-transformed Intensity", normType)
      ##   main <- sprintf("%s-transformed Mean/Bayes-regularized Dependence", normType)
      ## }
      ## print(xyplot(bayessd ~ normmean| type, toPlot, ylab='Bayes-regularized Standard Deviation',
      ##              xlab=xlab,
      ##              main=main))
      ## graphics.off()
      if (!any(is.na(toPlot$normsd))){
        cat('About to print bayessdmean_norm_plot..!\n')
        cat(sprintf('toPlot Head : %s', head(toPlot)), '\n')
        cat(sprintf('toPlot Dim : %s', dim(toPlot)), '\n')
        cat(sprintf('toPlot colnames : %s', colnames(toPlot)), '\n')
        graphicsFunc(sprintf('bayesvempsd_norm_plot.%s', fileNameFmt))
        main <- mainFunc("Raw Data Regularization Effect")
        if (normType != 'None'){
          main <- mainFunc(sprintf("%s-transformed Data Regularization Effect", normType))
        }
        print(xyplot(bayessd ~ normsd |type, toPlot, xlab='Empirical Standard Deviation',
                     ylab='Bayes-regularization Standard Deviation',
                     main=main,
                     panel = function(...) {
                       getPanelFunc()(...);
                       panel.abline(a = c(0, 1), lty = 2);
                     }, subset=!outlierFunc(toPlot$bayessd, toPlot$normsd)))
        graphics.off()
      
        ## Combine the above two:
        typeVect <- c(levels(toPlot$type)[toPlot$type], levels(toPlot$type)[toPlot$type])
        ## And then add in a plot where we look at the mean vs std for bayes/non bayes
        nToPlot <- data.frame(mean=c(toPlot$normmean, toPlot$normmean),
                              sd=c(toPlot$normsd, toPlot$bayessd),
                              sdtype=rep(c('Empirical', 'Bayes-regularized'), each=nrow(toPlot)),
                              type=typeVect)
        xlab <- "Mean Raw Intensity"
        main <- mainFunc("Effect of regularization - Mean/Variance Dependence")
        if (normType != 'None'){
          xlab <- sprintf("Mean %s-transformed Intensity", normType)
        }
        cat('About to bayesandempsd_vs_mean_plot.!\n')
        graphicsFunc(sprintf('bayesandempsd_vs_mean_plot.%s', fileNameFmt))
        cat(sprintf('getPanelFunc returns : %s \n', head(panelFunc)) )
        print(xyplot(sd~ mean | sdtype + type, nToPlot, xlab=xlab,
                     ylab="Standard Deviation", main=main,
                     panel=panelFunc, subset=!outlierFunc(nToPlot$sd, nToPlot$mean)))
        graphics.off()
      }
    }
  }
  
  if (! is.null(resData$ROC.x) & ! is.null(resData$ROC.y)){
    rocData <- data.frame(ROC.x=resData$ROC.x, ROC.y=resData$ROC.y)
    cat('About to print roc_plot..!\n')
    graphicsFunc(sprintf('roc_plot.%s', fileNameFmt))
    print(xyplot(ROC.y ~ ROC.x, rocData, xlab = "False Positive Rate [FP/(FP+TN)]",
         ylab = "True Positive Rate [TP/(TP+FN)]",
         main = mainFunc("Receiver Operating Characteristic (ROC) Curve")))
    graphics.off()
    write.table(rocData, file='ROC.txt', sep='\t', col.names=T, row.names=F)
  }
  graphics.off()
}

################################################################################
#### One-Way ANOVA Code #####
## runAllBayesAnova - wrapper to do everything for the ANOVA analysis
################################################################################
runAllBayesAnova <- function(aData, numVec, bayes=1, winSize=101, conf=10, ppde=1, betaFit=1,
                     doMulttest=FALSE, bayesIntAnova=FALSE, doPostHoc=FALSE, postHocType="T"){
  ## Method to run all parts of the ANOVA analysis;
  ## From the actual bayesAnova calls to doing ppde, mult hyp testing, and post-hoc tests
  aovRes <- bayesAnova(aData, numVec, bayes, winSize, conf, bayesIntAnova);
  if (doPostHoc){
    postHocData <- postHoc(aovRes, numVec, postHocType, bayes)
    aovRes <- data.frame(aovRes, postHocData)
  }
  if (ppde){
    ppdeData <- ppdeMix(as.matrix(aovRes$pVal), betaFit)
    aovRes <- data.frame(aovRes, ppdeData)
  }
  if (doMulttest){
    mtData <- runMulttest(aovRes$pVal)
    aovRes <- data.frame(aovRes, mtData)
  }
  aovRes
}

################################################################################
## bayesAnova
##
## adata is organized as: cond1-1,cond1-2,...,condN-1,condN-2,condN-m
## numVec is the number of replicates in each condition e.g. c(3,3,4)
## the number of columns should be the sum of numVec
## ppde, betaFit - ppde parameters
## doMulttest - do multiple testing
## bayesIntAnova - pool the intensities rather than the variances.
################################################################################
bayesAnova<-function(aData,numVec,bayes=1,winSize=101,conf=10, bayesIntAnova=F){
  numGenes <- dim(aData)[1]
  totMeasures <- dim(aData)[2]
  
  if (!(totMeasures==sum(numVec)))
    stop("ERROR: Number of replicates in each condition do not match total number of columns in dataset.")
  if (ceiling((winSize - 1) / 2) != ((winSize - 1) / 2))
    stop("ERROR: winSize must be an odd number.")
  winSize = (winSize - 1) / 2
  
  numCond=dim(as.array(numVec))
  
  ## store nVec, means, standard deviations in arrays
  nVec=matrix(data=NA,nrow=numGenes,ncol=numCond)
  mVec=matrix(data=NA,nrow=numGenes,ncol=numCond)
  sdVec=matrix(data=NA,nrow=numGenes,ncol=numCond)
  cat('numCond: ', numCond, '\n')
  for (i in 1:numCond){
    ##compute number of valid entries for each gene
    cat('The expr : ', (1+sum(numVec[0:(i-1)])):sum(numVec[1:i]), '\n')
    nVec[,i] <- apply(as.matrix(aData[, (1+sum(numVec[0:(i-1)])):sum(numVec[1:i])]),
                      1, function(x) sum(!is.na(x)))
    
    ## compute means for valid entries
    mVec[,i] <- apply(as.matrix(aData[, (1+sum(numVec[0:(i-1)])):sum(numVec[1:i])]),
                      1, function(x) {
                        if (sum(!is.na(x)))
                          mean(x[!is.na(x)])
                        else NA}
                      )
    ## check that at least 2 replicates are present
    sdVec[,i] <- apply(as.matrix(aData[, (1+sum(numVec[0:(i-1)])):sum(numVec[1:i])]),
                       1, function(x) {
                         if (sum(!is.na(x)) > 1)
                           sqrt(var(x[!is.na(x)]))
                         else NA}
                       )
    colnames(nVec)=paste("Num",1:numCond,sep="")
    colnames(mVec)=paste("Mean",1:numCond,sep="")
    colnames(sdVec)=paste("SD",1:numCond,sep="")
  }

  ## Calculate the statistics either in the bayesian or classical way
  if (bayes){
    ## If there is one value or explicitly told so, we should use the intensities
    if (any(numVec == 1)){bayesIntAnova = TRUE}
    
    rasdVec=matrix(data=NA,nrow=numGenes,ncol=numCond)
    colnames(rasdVec)=paste("bayesSD",1:numCond,sep="")

    if (!bayesIntAnova){
      ## compute rasd for each condition - using just averages
      ## over sd windows
      for (j in 1:numCond){
        stdC<-sdVec[,j]
        meanC<-mVec[,j]
        nC<-nVec[,j]
        temp <- runavg(stdC[!is.na(stdC)][order(meanC[!is.na(stdC)])],winSize)	
        temp <- temp[rank(meanC[!is.na(stdC)])]	
        rasdC <- rep(NA, numGenes)
        rasdC[!is.na(stdC)]<- temp
        
        ##compute bayes sd	
        rasdVec[,j]<- sqrt((conf * rasdC^2 + (nC - 1) * stdC^2)/(conf + nC - 2))
      }
    } else {
      ## Compute rasd for each condition
      ## Using sd over windows of intensities
      for (j in 1:numCond){
        stdC<-sdVec[,j]
        meanC<-mVec[,j]
        nC<-nVec[,j]
        intMat <- as.matrix(aData[!is.na(meanC), (1+sum(numVec[0:(j-1)])):sum(numVec[1:j])])
        intMat <- intMat[order(meanC[!is.na(meanC)]),]
        temp <- runavgPool(intMat, winSize)	
        temp <- temp[rank(meanC[!is.na(meanC)])]	
        rasdC <- rep(NA, numGenes)
        rasdC[!is.na(meanC)]<- temp

        ## Then make sure that meanC, stdC, and rasdC are all
        ## almost zero (if na)
        stdC[is.na(stdC) & !is.na(rasdC)] <- ALMOST_ZERO
        rasdC[rasdC == 0] <- ALMOST_ZERO
        
        ##compute bayes sd	
        rasdVec[,j]<- sqrt((conf * rasdC^2 + (nC - 1) * stdC^2)/(conf + nC - 2))
      }
    }

    aovBayes <- matrix(data=NA, nrow=numGenes, ncol=6)
    colnames(aovBayes) <- c("MSE.B", "MSE.W", "Fstat", "dfBetBayes", "dfWithBayes", "pVal")
    for (i in 1:numGenes){
      aovBayes[i,] <- bayesFtest(nVec[i,],
                                 mVec[i,],
                                 rasdVec[i,],
                                 conf=conf,bayes=1)
      if (!is.na(aovBayes[i,3])){
        if (aovBayes[i,3] < 0) aovBayes[i,1:6] <- NA 
      }
      else aovBayes[i,1:6] <- NA
    }
    aovBayes <- cbind(aData, nVec, mVec, sdVec, rasdVec, aovBayes)
    rm(nVec,mVec,sdVec,rasdVec)
  }
  else{
    ## Non-bayesian way of calculating
    aovBayes <- matrix(data=NA, nrow=numGenes, ncol=6)
    colnames(aovBayes) <- c("MSE.B", "MSE.W", "Fstat", "dfBet", "dfWith", "pVal")
    for (i in 1:numGenes){
      ## check that at least 2 replicates are present
      ## This must be true for the non-bayesian
      aovBayes[i,] <- bayesFtest(nVec[i,(nVec[i,]>1)],
                                 mVec[i,(nVec[i,]>1)],
                                 sdVec[i,(nVec[i,]>1)],
                                 conf=conf,bayes=0)
      if (!is.na(aovBayes[i,3])){
        if (aovBayes[i,3]<0) aovBayes[i,1:6] <- NA
      }
      else aovBayes[i,1:6] <- NA
    }
    aovBayes <- cbind(aData, nVec, mVec, sdVec, aovBayes)
    rm(nVec, mVec, sdVec)
  }
  return(aovBayes)
}


################################################################################
## Helper function for the ANOVA analysis
##
## Note: To properly handle the one sample case, we need to
## pass in a flag for the bayes
################################################################################
postHoc<-function(aovData, numVec, postTest="T", bayes=FALSE){
  ##postTest: postHoc test to be used, "T"=TukeyHSD, "S"=Scheffe's	
  numGenes <- dim(aovData)[1]
  numCond <- dim(as.array(numVec))
  lastCol <- dim(aovData)[2]
  meanCol <- (sum(numVec)+numCond+1):(sum(numVec)+2*numCond)
  numCol <- (sum(numVec)+1):(sum(numVec)+numCond)
  
  ##MSE (Within)
  mseCol <- lastCol-4
  dfBetCol <- lastCol-2
  dfWithCol <- lastCol-1
  
  postHocRes <- matrix(data=NA,nrow=nrow(aovData),ncol=(numCond*(numCond-1)/2))
  print(paste("Note: This analysis may be slow. Performing post-hoc tests for",
              numGenes,"genes in each of the",numCond*(numCond-1)/2,
              "pairwise comparisons.",sep=" "))
  colSt <- NULL	
  comparisonNum <- 0
  for (i in 1:numCond)
    for (j in i:numCond)
      if (i!=j){
        print(paste("Performing pairwise comparison",i,"-",j,sep=" "))
        comparisonNum <- comparisonNum+1
        colSt <- c(colSt,paste("Pair",i,"_",j,sep=""))
        if (postTest=="T")  
          postHocRes[,comparisonNum] <- tukeyhsd(aovData[,numCol[i]], aovData[,numCol[j]],
                                                 aovData[,meanCol[i]], aovData[,meanCol[j]],
                                                 aovData[,mseCol], aovData[,dfBetCol],
                                                 aovData[,dfWithCol], bayes)
        else if (postTest=="S")
          postHocRes[,comparisonNum] <- scheffes(aovData[,numCol[i]],aovData[,numCol[j]],
                                                 aovData[,meanCol[i]], aovData[,meanCol[j]],
                                                 aovData[,mseCol], aovData[,dfBetCol],
                                                 aovData[,dfWithCol], bayes)
      }
  rownames(postHocRes) <- rownames(aovData)
  colnames(postHocRes) <- colSt
  return(postHocRes)
}


scheffes <- function(nA, nB, mA, mB, MSE, dfBet, dfWith, bayes=FALSE){
  ## If bayes and number of replicates are one.., increase by one.
  if (bayes && nA == 1){nA <- 2}
  if (bayes && nB == 1){nB <- 2}
  S <- NULL
  S <- ((mA - mB)^2) / (MSE * (1 / nA + 1 / nB) * dfBet)
  ##check that at least 2 replicates are present				
  S[nA<=1] <- NA
  S[nB<=1] <- NA
  return(1 - pf(S, dfBet, dfWith))
}

tukeyhsd <- function(nA, nB, mA, mB, MSE, dfBet, dfWith, bayes=FALSE){
  if (bayes && nA == 1){nA <- 2}
  if (bayes && nB == 1){nB <- 2}
  ## Warning: When group sizes are unequal, the harmonic mean is used	
  hsd <- NULL
  hsd <- abs(mA - mB) / sqrt(MSE / (2 / (1 / nA + 1 / nB)))
  ##check that at least 2 replicates are present				
  hsd[nA<=1] <- NA
  hsd[nB<=1] <- NA
  return(1-ptukey(hsd,dfBet+1,dfWith))
}


bayesFtest<-function(nVec,mVec,sdVec,conf,bayes=1){
  ## nVec are numbers of control and expt
  ## mVec are the means resp.
  ## sdVec are the standard deviations resp.

  ## Note: if bayes == 1 AND nVec[i] == 1, we need to tweak
  ## the counts in the SS calcs.
  numGps <- dim(as.array(nVec))

  ## Calc the total group mean.
  gdMean <- 0
  for (i in 1:numGps)
    gdMean <- gdMean + nVec[i] * mVec[i]
  gdMean <- gdMean / sum(nVec)

  ## Calculate the explained (between-group) variance
  betSS <- 0
  for (i in 1:numGps)
    betSS <- betSS + nVec[i] * (mVec[i] - gdMean)^2

  ## Note the tweak we put in here if the nVec is 1 (and we are
  ## doing the bayesian analysis)
  withSS <- 0
  for (i in 1:numGps)
    withSS <- withSS + sdVec[i]^2 * (nVec[i] - 1 + (bayes && nVec[i] == 1))
  
  ## Regular ANOVA
  dfBet <- numGps - 1  			 			# k-1
  dfWith <- sum(nVec) - numGps					# k(n-1), n is mean nVec

  ## This will only happen if all singletons in the ANOVA.
  if (dfWith < numGps){
    dfWith <- numGps
  }
  
  dfTot <- sum(nVec) - 1					# kn-1
  ## f-statistic is computing on original dfs similar to t-test is computed on original n's
  fstat <- round((betSS / dfBet) / (withSS / dfWith), 4)
  pVal <- 1 - pf(fstat, dfBet, dfWith)
  
  ## Bayes ANOVA (Augmented degres of freedom, taking into account "conf")
  dfBetBayes <- numGps - 1           			        # k-1
  dfWithBayes <- sum(nVec) + numGps * conf - 2 * numGps		# k(n-1+conf-1), n is mean nVec
  dfTotBayes <- sum(nVec) + numGps * conf - numGps - 1		# k(n+(conf-1))-1
  pValBayes <- 1 - pf(fstat, dfBetBayes, dfWithBayes)
  
  if (bayes){
    return(cbind(betSS/dfBet, withSS/dfWith, fstat, dfBetBayes, dfWithBayes, pValBayes))
  }
  else{
    return(cbind(betSS/dfBet, withSS/dfWith, fstat, dfBet, dfWith, pVal))
  }
}


################################################################################
## Code for ppde calculation
################################################################################
ppdeMix <- function(pVal,n){
  ##Implements the method described in Allison et. al. 
  ##Computational Statistics & Data Analysis, 39:1-20 (2002) 
  
  ## The following objects are suggested objects for input into the mixture fitting
  ## They are lower and upper parameter bounds for the mixture parameters for a one beta,
  ## two beta, and three beta mixture, respectively.  The order of the parameter bounds is
  ## always, lambda0, lambda1, etc, then r1, s1, r2, s2, etc., depending on whether there
  ## are 1, 2, or 3 beta distributions.
  ## Note: beta parameters are capped at 170 since this value in the distribution,
  ## gamma(170) results in Splus reaching infinity limits.
  
  low1 <- c(.001,.01,.01)
  up1 <- c(1,165,165)
  low2 <- c(.001,.001,.01,.01,.01,.01)
  up2 <- c(1,1,165,165,165,165)
  low3 <- c(.001,.001,.001,.01,.01,.01,.01,.01,.01)
  up3 <- c(1,1,1,170,170,170,170,170,170)
  
  ## The following are suggested starting values for the mixture fitting
  ## The starting values are for parameters in the same order as above

  p01 <- c(.6,1,4)
  p02 <- c(.9,.5,4,1,1,4)
  p03 <- c(.8,.06,.1,.4,40,.4,55,1,6)
  
  if (ncol(pVal)==2){
    Pdata <- pVal[,2]
    Pdata <- Pdata[!is.na(Pdata)]
  }
  else
    Pdata <- pVal[!is.na(pVal)]
  Pdata<-adjPVal(Pdata)
  
  postP1 <- rep(NA, nrow(pVal))
  postP2 <- rep(NA, nrow(pVal))
  if(n==1){
    parms <- mle.mix(p01,mix.obj1,low1,up1,Pdata)
    L0 <- parms[1]
    L1 <- 1 - parms[1]
    r1 <- parms[2]
    s1 <- parms[3]
    
    numer1 <- L1 * pbeta(Pdata,r1,s1)
    denom1 <- numer1 + L0 * Pdata
    numer2 <- L1 * dbeta(Pdata,r1,s1)
    denom2 <- numer2 + L0
    
    rocX <- rep(NA, nrow(pVal))
    rocY <- rep(NA, nrow(pVal))
    if (ncol(pVal)==2){
      rocX[!is.na(pVal[,2])]<-Pdata
      rocY[!is.na(pVal[,2])]<-pbeta(Pdata,r1,s1)
    }
    else{
      rocX[!is.na(pVal)]<-Pdata
      rocY[!is.na(pVal)]<-pbeta(Pdata,r1,s1)
    }
  }
  else if(n==2){
    parms <- mle.mix(p02,mix.obj2,low2,up2,Pdata)
    L0 <- parms[1]
    L1 <- parms[2]*(1-parms[1])
    L2 <-(1-parms[1])*(1-parms[2])
    r1 <- parms[3]
    s1 <- parms[4]
    r2 <- parms[5]
    s2 <- parms[6]
    numer1 <- L1 * pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2)
    denom1 <- numer1 + L0 * Pdata
    numer2 <- L1 * dbeta(Pdata,r1,s1) + L2 * dbeta(Pdata,r2,s2)
    denom2 <- numer2 + L0
    
    rocX <- rep(NA, nrow(pVal))
    rocY <- rep(NA, nrow(pVal))
    if (ncol(pVal)==2){
      rocX[!is.na(pVal[,2])]<-Pdata
      rocY[!is.na(pVal[,2])]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2))/(L1+L2)
    }
    else{
      rocX[!is.na(pVal)]<-Pdata
      rocY[!is.na(pVal)]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2))/(L1+L2)
    }
  }
  else if(n==3){
    parms <- mle.mix(p03,mix.obj3,low3,up3,Pdata)
    L0 <- parms[1]
    L1 <- parms[2]*(1-parms[1])
    L2 <- parms[3]*(1-parms[1])*(1-parms[2])
    L3 <- (1-parms[3])*(1-parms[1])*(1-parms[2])
    r1 <- parms[4]
    s1 <- parms[5]
    r2 <- parms[6]
    s2 <- parms[7]
    r3 <- parms[8]
    s3 <- parms[9]
    numer1 <- L1 * pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2) + L3 * pbeta(Pdata,r3,s3)
    denom1 <- numer1 + L0 * Pdata
    numer2 <- L1 * dbeta(Pdata,r1,s1) + L2 * dbeta(Pdata,r2,s2) + L3 * dbeta(Pdata,r3,s3)
    denom2 <- numer2 + L0
    
    rocX <- rep(NA, nrow(pVal))
    rocY <- rep(NA, nrow(pVal))
    if (ncol(pVal)==2){
      rocX[!is.na(pVal[,2])]<-Pdata
      rocY[!is.na(pVal[,2])]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2) + L3 * 
                               pbeta(Pdata,r3,s3))/(L1+L2+L3)
    }
    else{
      rocX[!is.na(pVal)]<-Pdata
      rocY[!is.na(pVal)]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2) + L3 * 
                           pbeta(Pdata,r3,s3))/(L1+L2+L3)
    }
  }
  else cat("Error: Must fit 1,2, or 3 beta's (change betaFit)\n")
  if (ncol(pVal)==2){
    postP1[!is.na(pVal[,2])] <- numer1/denom1
    postP2[!is.na(pVal[,2])] <- numer2/denom2
  }
  else{
    postP1[!is.na(pVal)] <- numer1/denom1
    postP2[!is.na(pVal)] <- numer2/denom2
  }
  ##if (n==1)
  postP<- cbind(postP1,postP2,rocX,rocY)
  colnames(postP) <- c("cum.ppde.p","ppde.p","ROC.x","ROC.y")
  ##else
  ##	postP<- cbind(postP1,postP2)
  if (ncol(pVal)==2)
    postP <- cbind(pVal, postP)
  return(postP)
}

mix.obj1 <- function(p,x){
  ## This object is the objective function for a mixture of a uniform
  ## and one beta distribution
  e <- p[1] + (1-p[1])*(gamma(p[2]+p[3])/(gamma(p[2])*gamma(p[3])))*x^(p[2]-1)*(1-x)^(p[3]-1)
  sle <- -sum(log(e))
  sle
}

mix.obj2 <- function(p,x){
  ## This object is the objective function for a mixture of a uniform
  ## and two beta distributions
  e <- p[1] + (p[2]*(1-p[1]))*(gamma(p[3]+p[4])/(gamma(p[3])*gamma(p[4])))*x^(p[3]-1)*(1-x)^(p[4]-1) +
    ((1-p[1])*(1 - p[2]))*(gamma(p[5]+p[6])/(gamma(p[5])*gamma(p[6])))*x^(p[5]-1)*(1-x)^(p[6]-1)
  sle <- -sum(log(e))
  sle
}

mix.obj3 <- function(p,x){
  ## This object is the objective function for a mixture of a uniform
  ## and three beta distributions
  e <- p[1] + (p[2]*(1-p[1]))*(gamma(p[4]+p[5])/(gamma(p[4])*gamma(p[5])))*x^(p[4]-1)*(1-x)^(p[5]-1) +
    (p[3]*(1-p[1])*(1 - p[2]))*(gamma(p[6]+p[7])/(gamma(p[6])*gamma(p[7])))*x^(p[6]-1)*(1-x)^(p[7]-1) +
      ((1-p[3])*(1-p[1])*(1 - p[2]))*(gamma(p[8]+p[9])/(gamma(p[8])*gamma(p[9])))*x^(p[8]-1)*(1-x)^(p[9]-1)
  sle <- -sum(log(e))
  sle
}

mle.mix <- function(init,mix.obj,low,up,Pdata){
  ## This object computes the MLE's for the mixture distribution
  ## Inputs are:
  ## init = starting values for the computations (see p01, p02, p03 above)
  ## mix.obj = is the objective function (i.e., see mix.obj1 etc from above)
  ## low and up are the parameter limits (see low1, up1 etc from above)
  ## Pdata = the p-values.  This should be a single vector of length
  ## equal to the number of genes in the study.
  ## Note: all inputs should be consistent with regards to the number of components
  ## being fitted.  That is, p01 goes with low1, up1, and mix.obj1.
  ## Output:  An object called mix that is the result of the fitting algorithm.
  ## mix$par are the MLE's of parameters in the same order as above objects.
  ## Note: the weighting parameters on mixture components need to be reparametrized
  ## to make sense.
  tryCatch({ mix <- optim(init,mix.obj,lower=low, upper=up, method="L-BFGS-B", x=Pdata);},
           error = function(e){
             cat('****Fatal CyberT ERROR****\n',
                 'Unable to fit the PPDE Mixture Model.\n',
                 'Try not running the PPDE or increase the confidence/window size.\n',
                 '****End Fatal CyberT ERROR****\n')
             stop(e)}
             )
  return(mix$par)
}

## Some PPDE Functions for fitting and analysis:
## Mainly these are some ways to look at different ways of fitting
## 2 beta + 1 uniform distributions.
## 4-14-2012 - Not used in production code, but available to be used
##  from R directly, or incorporated into site at a later date.

plotBetaFromParm <- function(parm){
  ## Given a param vector (as from mle.mix) for 2 beta setup,
  ## plot the pdfs of the two fit betas.
  library(lattice)
  r1 <- parm[3]; s1 <- parm[4]
  r2 <- parm[5]; s2 <- parm[6]
  
  x <- seq(0.001,.999,.01)
  toplot <- data.frame(x=x, b1=dbeta(x, r1, s1), b2=dbeta(x, r2, s2) )
  xyplot(b1 + b2 ~ x, toplot, type='l',
         par.setting=list(superpose.line=list(col=c('red', 'blue'))),
         auto.key=list(points=F, lines=T, corner=c(0.2, .6), size=3))
} 

dgamma <- function(x){
  ## Derivative of gamma function at x:
  digamma(x)*gamma(x)
}

## dbeta.a, dbeta.b - partial derivs of beta function for a,b
dbeta.a <- function(a, b){
  (gamma(a + b)*gamma(b)*dgamma(a) - gamma(a)*gamma(b)*dgamma(a + b))/(gamma(a + b))^2
}
dbeta.b <- function(a, b){
  (gamma(a + b)*dgamma(b)*gamma(a) - gamma(a)*gamma(b)*dgamma(a + b))/(gamma(a + b))^2
}

em.mle.mix.2beta <- function(init,low,up,pdata){
  ## Alternative way to fit the two beta scheme using Expectation Maximization.
  ##
  
  piparams <- init[1:2]
  ## 1 - (sum(piparams)) <- second beta weight

  thetaparams <- init[-c(1:2)]
  theta.low <- low[-c(1:2)]
  theta.up <- up[-c(1:2)]
  
  ## Now make a matrix of weights for each data point
  weights <- matrix(rep(c(piparams[1], piparams[2], (1- piparams[1] - piparams[2])),
                        each=length(pdata)), nrow=length(pdata))
  
  updateWeights <- function(w, thetas){
    ## E-Step - given some thetas, update the weights.
    ## This is analytic, so simple
    w[, 1] <- 1
    w[, 2] <- dbeta(pdata, thetas[1], thetas[2])
    w[, 3] <- dbeta(pdata, thetas[3], thetas[4])
    sw <- apply(w, 1, sum)
    w <- w/sw
    w
  }

  objFunc <- function(thetas, w){
    ## Expected complete-data log-likelihood (to optimize in M-step)
    ## Calc the beta vals, but give a little weight to 
    b1 <- dbeta(pdata, thetas[1], thetas[2])
    b1[b1 < .Machine$double.eps] <- .Machine$double.eps
    b2 <- dbeta(pdata, thetas[3], thetas[4])
    b2[b2 < .Machine$double.eps] <- .Machine$double.eps
    
    r <- w[,2] * log(b1) +
      w[, 3] * log(b2)
    -sum(r)
  }
  d.objFunc <- function(p, w){
    ## Can break this down per theta    
    a1 <- w[, 2] * (-dbeta.a(p[1], p[2])/beta(p[1], p[2]) +
                    log(pdata))
    b1 <- w[, 2] * (-dbeta.b(p[1], p[2])/beta(p[1], p[2]) +
                    log(1-pdata))
    a2 <- w[, 3] * (-dbeta.a(p[3], p[4])/beta(p[3], p[4]) +
                    log(pdata))
    b2 <- w[, 3] * (-dbeta.b(p[3], p[4])/beta(p[3], p[4]) +
                    log(1-pdata))
    a1 <- -sum(a1)
    b1 <- -sum(b1)
    a2 <- -sum(a2)
    b2 <- -sum(b2)

    ## Getting some strange infinite sized steps, cut down
    ##  the gradient calc by a factor of 10
    0.1 * c(a1, b1, a2, b2)
  }
  updateThetas <- function(thetas, weights){
    ## M-Step - given some weights, update the thetas.
    ## This cannot be solved analytically, so have to do some optimization.
    opt <- optim(thetas, objFunc, d.objFunc, lower=theta.low, upper=theta.up,
                 method="L-BFGS-B", w=weights)
    opt$par
  }

  ## Initial E-Step
  weights <- updateWeights(weights, thetaparams)

  minDiff <- 1
  niter <- 0
  while (minDiff > 2*sqrt(.Machine$double.eps)){
    old.thetas <- thetaparams
    old.objVal <- objFunc(thetaparams, weights)
    ## M-Step
    thetaparams <- updateThetas(thetaparams, weights)

    ## E-Step 
    weights <- updateWeights(weights, thetaparams)

    ## for debugging
    niter <- niter + 1
    objVal <- objFunc(thetaparams, weights)
    #cat('At Iter:', niter, ' with objFunc: ', objVal, '\n')
    
    ## calc divergence
    minDiff <- abs(old.objVal - objVal)
    if (niter > 100){
      minDiff <- .Machine$double.eps
    }
  }

  final.weights <- updateWeights(weights, thetaparams)
  pivals <- apply(final.weights, 2, sum)/nrow(final.weights)

  return( c(pivals[1], pivals[2], thetaparams))
}

################################################################################
## End of ppde section
################################################################################

writeError <- function(msg){
  cat('****Fatal CyberT ERROR****\n',
        msg, '\n',
        '****End Fatal CyberT ERROR****\n')
}


