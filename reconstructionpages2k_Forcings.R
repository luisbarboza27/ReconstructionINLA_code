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


###OBSERVED TEMPERATURES
tempnh <- read.table(file = './data/HadCRUT.4.4.0.0.ama_ns_avg_1850-2015.txt')
names(tempnh)=c('year','temp')
tempobs <- tempnh %>% filter(year %in% seq(1900,2014)) 
Tnh <- tempobs %>% filter(year %in% seq(1900,2000)) %>% dplyr::select(temp)

###REDUCED PROXIES###########
load(paste0('./results/RPs/RP_new_',dataset,'_',method,'.RData'))
Rp <- RP[,paste0('RP',1)]

###FORCINGS#################
forcing <- read.csv(file = './Forcings/forcing.csv')

Ss <- forcing %>% dplyr::select(solar) %>% mutate(S=(solar-mean(solar))/sd(solar)) %>% 
  dplyr::select(S)
Vs <- forcing %>% dplyr::select(volcanic) %>% mutate(V = log(-volcanic+1)) %>% 
  mutate(V = (V-mean(V))/sd(V)) %>% dplyr::select(V)
Cs <- forcing %>% dplyr::select(CO2) %>% mutate(C = log(CO2)) %>% 
  mutate(C = (C-mean(C))/sd(C)) %>% dplyr::select(C)


#####INLA FITTING
Y <- matrix(NA, N+N, 2)
Y[1:n,     1] <- Rp              # actual observations
Y[1:N + N, 2] <- 0              # faked observations

prior.beta <- c(0, 3)
prior.alpha0 <- c(0, 3)
prior.prec.epsilon <- c(1,1e-20)
prior.prec.eta <- c(1,1e-20)
inla.data = list() 

inla.data$Y <- Y
inla.data$Ss <- c(rep(NA,N),Ss$S,rep(0,mfore)) 
inla.data$Vs <- c(rep(NA,N),Vs$V,rep(0,mfore))
inla.data$Cs <- c(rep(NA,N),Cs$C,rep(5.9,mfore))
inla.data$alpha0 <- c(rep(1,N),rep(NA,N)) 
inla.data$beta0 <- c(rep(NA,N),rep(1,N)) 
inla.data$Ttime <- c(rep(NA,N),seq(1,N))
inla.data$Tweights <- c(rep(NA,N),rep(-1,N))
inla.data$alpha1time <- c(seq(1,N),rep(NA,N))
inla.data$etatime <- c(rep(NA,N), 1:N) 
inla.data$etaweights <- c(rep(NA,N), rep(1,N)) 

formula <- Y ~ -1+alpha0+beta0+Ss+Vs+Cs+
  f(alpha1time,copy = 'Ttime',hyper = list(beta = list(param = prior.beta, fixed = FALSE)))+
  f(Ttime,Tweights,model = 'iid',fixed = T,initial = -60,constr = F)+
  f(etatime,etaweights,model = 'iid',hyper=list(theta=list(param=prior.prec.eta)))


r <- inla(formula,data = inla.data,family = c('gaussian','gaussian'),
          control.family = list(
            list(hyper=list(theta = list(initial = 10,
                                         param = prior.prec.epsilon,
                                         fixed = FALSE))),
            list(hyper=list(theta=list(initial = 60,fixed = TRUE)))),
          control.fixed = list(
            mean = list(alpha0=prior.alpha0[1],beta0=prior.beta[1],Ss=prior.beta[1],Vs=prior.beta[1],Cs=prior.beta[1]),
            prec = list(alpha0=prior.alpha0[2],beta0=prior.beta[2],Ss=prior.beta[2],Vs=prior.beta[2],Cs=prior.beta[2]))
)

data1899 <- data.frame(Tt = r$summary.random$Ttime$mean[1:1899],
                       Tt025 = r$summary.random$Ttime$`0.025quant`[1:1899],
                       Tt975 = r$summary.random$Ttime$`0.975quant`[1:1899])
data2000 <- data.frame(Tt = r$summary.random$Ttime$mean[1900:2000],
                       Tt025 = r$summary.random$Ttime$`0.025quant`[1900:2000],
                       Tt975 = r$summary.random$Ttime$`0.975quant`[1900:2000])


save.image(file = paste0('./results/results_new_',dataset,'_',method,'_',1,'_forcings.RData'))
