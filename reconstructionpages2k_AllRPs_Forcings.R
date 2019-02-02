##THIS SCRIPT COMPUTES THE PALEOCLIMATE RECONSTRUCTION ACCORDING TO
## BARBOZA ET AL, 2019. (SCENARIO WF)

library(INLA)
library(dplyr)
library(ggplot2)
library(stringr)

###PRELIMINARIES
horizonte <- seq(1,2000)
mfore <- 0
n <- length(horizonte)
N <- n+mfore
lagN <- 250
lenghts <- c(seq(0,nRPs-1),0)
Nvector <- N-lenghts*lagN


###OBSERVED TEMPERATURES
tempnh <- read.table(file = './data/HadCRUT.4.4.0.0.ama_ns_avg_1850-2015.txt')
names(tempnh)=c('year','temp')
tempobs <- tempnh %>% filter(year %in% seq(1900,2014)) 
Tnh <- tempobs %>% filter(year %in% seq(1900,2000)) %>% select(temp)

###REDUCED PROXIES###########
load(paste0('./results/RPs/RP_new_',dataset,'_',method,'.RData'))
Rp <- RP[,paste0('RP',1:8)]

###FORCINGS#################
forcing <- read.csv(file = './Forcings/forcing.csv')

Ss <- forcing %>% dplyr::select(solar) %>% mutate(S=(solar-mean(solar))/sd(solar)) %>% 
  dplyr::select(S)
Vs <- forcing %>% dplyr::select(volcanic) %>% mutate(V = log(-volcanic+1)) %>% 
  mutate(V = (V-mean(V))/sd(V)) %>% dplyr::select(V)
Cs <- forcing %>% dplyr::select(CO2) %>% mutate(C = log(CO2)) %>% 
  mutate(C = (C-mean(C))/sd(C)) %>% dplyr::select(C)


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
inla.data$Ss <- c(rep(NA,sum(Nvector[-(nRPs+1)])),Ss$S,rep(0,mfore)) 
inla.data$Vs <- c(rep(NA,sum(Nvector[-(nRPs+1)])),Vs$V,rep(0,mfore))
inla.data$Cs <- c(rep(NA,sum(Nvector[-(nRPs+1)])),Cs$C,rep(5.9,mfore))

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


inla.data$beta0 <- c(rep(NA,sum(Nvector[-(nRPs+1)])),rep(1,N)) 
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

formula <- Y ~ -1+alpha01+alpha02+alpha03+alpha04+alpha05+alpha06+alpha07+alpha08+beta0+Ss+Vs+Cs+
  f(alpha11time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha12time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha13time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha14time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha15time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha16time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha17time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(alpha18time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(Ttime,Tweights,model = 'iid',fixed = T,initial = -60,constr = F)+
  f(etatime,etaweights,model = 'iid',hyper=list(theta=list(param=prior.prec.eta)))

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
            mean = list(alpha01=prior.alpha0[1],alpha02=prior.alpha0[1],alpha03=prior.alpha0[1],alpha04=prior.alpha0[1],alpha05=prior.alpha0[1],alpha06=prior.alpha0[1],alpha07=prior.alpha0[1],alpha08=prior.alpha0[1],beta0=prior.beta[1],Ss=prior.beta[1],Vs=prior.beta[1],Cs=prior.beta[1]),
            prec = list(alpha01=prior.alpha0[2],alpha02=prior.alpha0[2],alpha03=prior.alpha0[2],alpha04=prior.alpha0[2],alpha05=prior.alpha0[2],alpha06=prior.alpha0[2],alpha07=prior.alpha0[2],alpha08=prior.alpha0[2],beta0=prior.beta[2],Ss=prior.beta[2],Vs=prior.beta[2],Cs=prior.beta[2]))
          ,verbose = T)



save.image(file = paste0('results_new_',dataset,'_',method,'_',nRPs,'_forcings.RData'))
