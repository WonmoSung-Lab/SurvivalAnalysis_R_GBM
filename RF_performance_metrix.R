# Import packages
library(randomForestSRC) # : random forest package
library(caret) # : to stratified k-fold sets
library(survival) # : to make survival object
library(mlr) # : to tune hyper parameters etc.
library(pec) # : to get IBS value
library(risksetROC)

# 1) Cindex
get_Cindex = function(pred_results, test_df){
  Cindex_all = c()
  for (time_i in 1:length(pred_results$time.interest)){
    #print(pred_results$survival[, time_i])
    Cindex_i = randomForestSRC::get.cindex(test_df$time, test_df$event, pred_results$survival[, time_i])
    Cindex_all = append(Cindex_all, Cindex_i)
  }
  Cindex_mean = mean(Cindex_all)
  return(Cindex_mean)
}

# 2) IBS
get_IBS <- function(rfsrc.obj, self=1, ave.flag=1){
  obj <- rfsrc.obj
  nsample <- obj$n
  time <- obj$time.interest
  ntime <- length(time)
  
  deadStatus <- obj$yvar[,2]
  deadTime <- obj$yvar[deadStatus==1,1]
  censorTime <- obj$yvar[deadStatus==0,1]
  
  if(self==0){
    survP <- obj$survival.oob
  }
  else{
    if(self==1){ # Not oob ibs. 
      survP <- obj$survival
    }}
  
  if(length(survP) == (nsample*ntime)){
    
    tmp <- matrix(survP, ncol=ntime)
    tmp1 <- tmp[deadStatus==1, ] # dead pt
    tmp0 <- tmp[deadStatus==0, ] # censored pt
    bs1 <- matrix(nrow=nrow(tmp1), ncol=ntime)
    bs0 <- matrix(nrow=nrow(tmp0), ncol=ntime)
    
    # go through dead 
    for(si in 1:nrow(tmp1)){ # per patient. 
      coreBefore <- c(1:ntime)[ time < deadTime[si] ]
      coreAfter <- c(1:ntime)[ time >= deadTime[si] ]
      bs1[si, coreBefore] <- (1-tmp1[si,coreBefore])^2
      bs1[si, coreAfter] <- tmp1[si,coreAfter]^2
    }
    
    # go through censored
    for(si in 1:nrow(tmp0)){
      coreBefore <- c(1:ntime)[ time < censorTime[si] ]
      coreAfter <- c(1:ntime)[ time >= censorTime[si] ]
      bs0[si, coreBefore] <- (1-tmp0[si,coreBefore])^2
      # it should not be zero, but unknown (not contributing)
      # bs0[si, coreAfter] <- 0
      bs0[si, coreAfter] <- NA
    }
    
    bs <- rbind(bs1, bs0)
    ibs.time1 <- apply(bs1, 2, mean)
    ibs.time0 <- apply(bs0, 2, function(x){ mean( sort(x))  })
    ibs.time <- apply(bs, 2, function(x){ mean( sort(x))  } )
    ibs.time.all <- rbind(ibs.time1,ibs.time0,ibs.time)
    
    if(ave.flag==1){
      aveibs <- mean(ibs.time.all[3,])
      return(aveibs)
    } else if(ave.flag==0){
      return(ibs.time.all)
    }
    
  }}




