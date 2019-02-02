##THIS SCRIPT COMPUTES THE PALEOCLIMATE RECONSTRUCTION ACCORDING TO
## BARBOZA ET AL, 2019. (SCENARIO NF)

library(INLA)
library(dplyr)
library(ggplot2)
library(fda)
library(stringr)



start_time <- Sys.time()
calibrationp <- 1900:2000

###PRELIMINARIES
horizonte <- seq(1,2000)
mfore <- 0
n <- length(horizonte)
N <- n+mfore
lagN <- 250
nRPs <- 8
lenghts <- c(seq(0,nRPs-1),0)
Nvector <- N-lenghts*lagN

nbasisv <- 6*(2000/100)
IS95 <- rep(0,length(nbasisv))
IS80 <- rep(0,length(nbasisv))
CRPSs <- rep(0,length(nbasisv))
sds <- rep(0,length(nbasisv))
bias <- rep(0,length(nbasisv))

###OBSERVED TEMPERATURES
tempnh <- read.table(file = './data/HadCRUT.4.4.0.0.ama_ns_avg_1850-2015.txt')
names(tempnh)=c('year','temp')
tempobs <- tempnh %>% filter(year %in% seq(1900,2014))
Tnh <- tempobs %>% filter(year %in% calibrationp) %>% select(temp)


###REDUCED PROXIES###########
load(paste0('./results/RPs/RP_new_',dataset,'_',method,'.RData'))
Rp <- RP[,paste0('RP',seq(1,8))]

for(j in 1:length(nbasisv)){
  show(j)
  BSplinesB <- create.bspline.basis(c(0,2001),nbasis = nbasisv[j])  #0-2001
  BSE <- eval.basis(evalarg = seq(1,2000),basisobj = BSplinesB)
  BSE <- as.data.frame(BSE)
  BSE <- scale(BSE)
  BSE <- cbind(rep(1,2000),BSE)
  nbasis <- nbasisv[j]+1
  
  
  #####INLA FITTING
  Y <- matrix(NA, sum(Nvector), nRPs+1)
  
  for(i in 1:nRPs){
    if(i ==1){
      Y[1:Nvector[i],i] <- Rp[,i]
    }else{
      Y[1:Nvector[i]+sum(Nvector[seq(1,i-1)]),i] <- Rp[!is.na(Rp[,i]),i]
    }
  }
  Y[1:N+sum(Nvector[seq(1,nRPs)]),nRPs+1] <- 0
  
  prior.beta <- c(0, 3)
  prior.alpha0 <- c(0, 3)
  prior.prec.epsilon <- c(1,1e-20)
  prior.prec.eta <- c(1,1e-20)
  inla.data = list()
  
  inla.data$Y <- Y
  nelements <- length(inla.data)
  for(i in 1:nbasis){
    inla.data[[nelements+i]] <- c(rep(NA,sum(Nvector[-(nRPs+1)])),BSE[,i],rep(0,mfore))
  }
  
  names(inla.data)[nelements+seq(1,nbasis)] <- paste0('BS',seq(1,dim(BSE)[2]))
  
  
  nelements <- length(inla.data)
  for(i in 1:nRPs){
    inla.data[[nelements+i]] <- rep(NA,sum(Nvector))
    if(i ==1){
      inla.data[[nelements+i]][1:Nvector[i]] <- rep(1,Nvector[i])
    }else{
      inla.data[[nelements+i]][1:Nvector[i]+sum(Nvector[seq(1,i-1)])] <- rep(1,Nvector[i])
    }
  }
  
  names(inla.data)[nelements+seq(1,nRPs)] <- paste0('alpha0',seq(1,nRPs))
  
  
  inla.data$Ttime <- c(rep(NA,sum(Nvector[-(nRPs+1)])),seq(1,N))
  inla.data$Tweights <- c(rep(NA,sum(Nvector[-(nRPs+1)])),rep(-1,N))
  
  nelements <- length(inla.data)
  for(i in 1:nRPs){
    inla.data[[nelements+i]] <- rep(NA,sum(Nvector))
    if(i ==1){
      inla.data[[nelements+i]][1:Nvector[i]] <- seq(1,Nvector[i])
    }else{
      inla.data[[nelements+i]][1:Nvector[i]+sum(Nvector[seq(1,i-1)])] <- seq((i-1)*lagN+1,N)
    }
  }
  
  names(inla.data)[nelements+seq(1,nRPs)] <- paste0('alpha1',seq(1,nRPs),'time')
  
  
  inla.data$etatime <- c(rep(NA,sum(Nvector[-(nRPs+1)])),1:N)
  inla.data$etaweights <- c(rep(NA,sum(Nvector[-(nRPs+1)])),rep(1,N))
  
  
  alphas0 <- str_c(paste0('alpha0',seq(1,nRPs)),collapse = '+')
  splines <- str_c(paste0('BS',seq(1,nbasis)),collapse = '+')
  alphas1 <- "f(alpha11time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha12time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha13time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha14time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha15time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha16time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha17time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha18time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(Ttime,Tweights,model = 'iid',fixed = T,initial = -60,constr = F)+
  f(etatime,etaweights,model = 'iid',hyper=list(theta=list(param=prior.prec.eta)))"
  
  rhs <- paste(alphas0,splines,alphas1,sep = '+')
  formula <- reformulate(termlabels = rhs,response = 'Y',intercept = F)
  
  
  fixed.prior.mean <- list()
  fixed.prior.sd <- list()
  
  for(i in 1:nRPs){
    fixed.prior.mean[[i]] <- prior.alpha0[1]
    fixed.prior.sd[[i]] <- prior.alpha0[2]
  }
  
  names(fixed.prior.mean)[seq(1,nRPs)] <- paste0('alpha0',seq(1,nRPs))
  names(fixed.prior.sd)[seq(1,nRPs)] <- paste0('alpha0',seq(1,nRPs))
  
  nelements <- length(fixed.prior.mean)
  
  for(i in 1:nbasis){
    fixed.prior.mean[[nelements+i]] <- prior.beta[1]
    fixed.prior.sd[[nelements+i]] <- prior.beta[2]
  }
  
  names(fixed.prior.mean)[nelements+seq(1,nbasis)] <- paste0('BS',seq(1,nbasis))
  names(fixed.prior.sd)[nelements+seq(1,nbasis)] <- paste0('BS',seq(1,nbasis))
  
  
  r <- inla(formula,data = inla.data,family = rep('gaussian',nRPs+1),
            control.family = list(
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta = list(initial = 10,
                                           param = prior.prec.epsilon,
                                           fixed = FALSE))),
              list(hyper=list(theta=list(initial = 60,fixed = TRUE)))),
            control.fixed = list(
              mean = fixed.prior.mean,
              prec = fixed.prior.sd)
            ,verbose = T)

}

save.image(file = paste0('./results/results_new_',dataset,'_',method,'_',nbasisv,'_Splines.RData'))
