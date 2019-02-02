##MASTER FILE THAT COMPUTES THE RECONSTRUCTION ACCORDING TO 
##BARBOZA ET AL, 2019. 

##Choose the scenario type (WF, NF, Mixed)
scenario='Mixed'

##Choose the reduced-proxy method (PCR, sPCR, SPLS, LASSO, SIR)
method='PCR' 

##Choose the dataset type (CW_FDR: screening-based dataset, All: without screening dataset)
dataset <- 'CW_FDR'

##Choose the number of reduced proxies to use. (Only 1 or 8)
nRPs <- 1


if(scenario=='WF'){
  if(nRPs==1){
    source('reconstructionpages2k_Forcings.R')
  }
  if(nRPs==8){
    source('reconstructionpages2k_AllRPs_Forcings.R')
  }
}

if(scenario=='NF'){
  source('reconstructionpages2k_AllRPs_BSplines.R')
}

if(scenario=='Mixed'){
  source('reconstructionpages2k_AllRPs_BSplines_mixed.R')
}

