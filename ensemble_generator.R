library(tidyverse)
library(INLA)

##Choose the number of samples (ensembles)
nsamples <- 10000

##Choose the time period to evaluate the samples
tperiod <- 1:2000

##Choose the scenario type (WF, NF, Mixed)
scenario='Mixed'

##Choose the reduced-proxy method (PCR, sPCR, SPLS, LASSO, SIR)
method='PCR' 

##Choose the dataset type (CW_FDR: screening-based dataset, All: without screening dataset)
dataset <- 'CW_FDR'

##Choose the number of reduced proxies to use. (Only 1 or 8)
nRPs <- 8

if(scenario=='WF'){
  load(file = paste0('./results/results_new_',dataset,'_',method,'_',nRPs,'_forcings.RData'))
}

if(scenario=='NF'){
  nbasisv <- 120
  load(file = paste0('./results/results_new_',dataset,'_',method,'_',nbasisv,'_Splines.RData'))
}

if(scenario=='Mixed'){
  nbasisv <- 100
  load(file = paste0('./results/results_new_',dataset,'_',method,'_',nbasisv,'_Splines_mixed.RData'))
}

samples_matrix <- map_dfc(tperiod,~inla.rmarginal(nsamples,marginal = r$marginals.random$Ttime[[.]]))
