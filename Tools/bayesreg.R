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
############################## bayesT ##################################################
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
  
  numGene<- base::nrow(aData)
  
  ## compute number of valid entries for each gene
  nC <- apply(base::as.matrix(aData[, 1:numC]), 1, function(x) sum(!is.na(x)))
  nE<- apply(base::as.matrix(aData[, (numC+1):(numC+numE)]), 1, function(x) sum(!is.na(x)))
  
  ## compute means for valid entries
  meanC<- apply(base::as.matrix(aData[, 1:numC]), 1, function(x) 
                if (sum(!is.na(x)))
                mean(x[!is.na(x)])
                else NA)
  meanE<- apply(base::as.matrix(aData[, (numC+1):(numC+numE)]), 1, function(x) 
                if (sum(!is.na(x)))
                mean(x[!is.na(x)])
                else NA)
  stdC<- apply(base::as.matrix(aData[, 1:numC]), 1, function(x) 
               if (sum(!is.na(x)) > 1)
               sqrt(var(x[!is.na(x)]))
               else NA)
  stdE<- apply(base::as.matrix(aData[, (numC+1):(numC+numE)]), 1, function(x) 
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
      intMat <- base::as.matrix(aData[!is.na(meanC), 1:numC]);
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
      intMat<- base::as.matrix(aData[!is.na(meanE), (numC+1):(numC+numE)]);
      intMat <- intMat[order(meanE[!is.na(meanE)]), ];
      temp <- runavgPool(intMat, winSize)	
      temp <- temp[rank(meanE[!is.na(meanE)])]
      rasdE[!is.na(meanE)]<- temp
    }
    
    ## Before computing the Bayes SD, go through and set StdC to almost (zero)  wherever is.na and rasdC is not NA
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
  pVal <- 1 - stats::pf(ttest[,1]^2, 1, ttest[,2])
  
  ## Then the final data.frame 
  if (bayes){
    objBayes<- cbind(pVal)
  }
  else{
    objBayes<- cbind(pVal)
  }

  ## Mult testing correction
  if (doMulttest){
    ##Bind them, last two cols of adjp list in original order
    objBayes <- data.frame(objBayes, runMulttest(pVal))
  }
  
  return(objBayes)
}

############################## runVsn #################################################
## Simple wrapper to load and then run vsn
################################################################################
runVsn <- function(data, ...){
  .libPaths("/home/baldig/shared_libraries/centos64/pkgs/R/2.15.1_fix/lib64/R/library")
  suppressMessages(require(vsn))
  #cat('data in VSN: \n', colnames(data), '\n', as.vector(unlist(data[1,])), '\n',
  #    as.vector(unlist(data[2,])), '\n', as.vector(unlist(data[3, ])), '\n');
  tryCatch(data.frame(vsn::justvsn(base::as.matrix(data), verbose=FALSE, minDataPointsPerStratum=10, ...)),
           error=function(e){
             libs_a<-.libPaths()
             writeError('Error in running vsn.  Possible collinearity.  Check your data.')
             stop(e)
           }
         )
}

############################## runMulttest ############################################
## Simple function to do Bonferroni and BH multiple testing correction.
################################################################################
runMulttest <- function(pvals){
  .libPaths("/home/baldig/shared_libraries/centos64/pkgs/R/2.15.1_fix/lib64/R/library")
  suppressMessages(require(multtest));
  adjPObj <- multtest::mt.rawp2adjp(pvals, proc=c("BH"))
  ##Bind them, last two cols of adjp list in original order
  adjPObj$adjp[order(adjPObj$index), -1]
}

#######################################################################################
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
  
  x <- base::as.matrix(x);
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

writeError <- function(msg){
  cat('****Fatal CyberT ERROR****\n',
        msg, '\n',
        '****End Fatal CyberT ERROR****\n')
}


