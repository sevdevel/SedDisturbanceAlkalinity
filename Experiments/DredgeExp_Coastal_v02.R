#---
#title: 'Impact chronic trawling on C and alkalinity generation'
#author: "Sebastiaan van de Velde"
#date: "12-06-2022"
#output:
#  pdf_document: default
#  html_document: default
#---

## Load the package for use
library(marelac)
library(viridis)

source('../Functions/TrawlImpactRecovery_script_v01.R')
source('../Functions/DiagCModel_v04.R')
source('../Functions/SolveSaphePackage.R')
source('../Functions/AuxFunctions_v04.R')

#======================================================================================================================
# Estimate annual dredging depth and frequency
#======================================================================================================================

# Global annual dredging activities may remove hundreds of millions to several billion cubic meters of sediment
# may involve dredging thousands to tens of thousands of hectares (1 hectare = 10,000 square meters) of surface area

V.dredge <- c(100.,1000.)*1e6   # [m3]
A.dredge <- c(1000.,10000.)*1e4 # [m2]
d.dredge <- V.dredge/A.dredge   # [m]
  
#
# load baseline file
load("../Output/Baseline_DepthTransect.Rdata")
PL     <- ModelResults[[1]]$PL
SV.ini <- ModelResults[[1]]$output$y

# Adapt disturbance parameters 
PL$yrs <- 60
# PD in mud= 3.2 +/- 1. 95%CI[1.4, 5.2]
# Average erosion Beam Trawl
PL$mud.content <- 50.
PL$depth.ero   <- 20.#sum((2.602*PL$mud.content/100+1.206e-3*c(13,1967,572,2118) + 1.321e-2*PL$mud.content/100*c(13,1967,572,2118))/(PL$rho.sed*1e3)*1e2) # Beam trawl large vessel
#PL$depth.ero <- sum((2.602*50/100+1.206e-3*c(19,1595,135,699) + 1.321e-2*50/100*c(19,1595,135,699))/(2.6*1e3)*1e2)  # Beam trawl small vessel
#PL$depth.ero <- sum((2.602*50/100+1.206e-3*c(1158,7.6,4020) + 1.321e-2*50/100*c(1158,7.6,4020))/(2.6*1e3)*1e2)      # otter trawl 'average Scottish' vessel
PL$depth.mix  <- 20.
PL$res.effect <- c(1.,1.,1.) 
PL$longevity   <- 7.5
PL$TI.Benthos  <- FALSE
# Update parameter list 
PL <- set.parameter.list(PL)

#create sensitivity range
sens.par <- list()
#sens.par$frequency <- c(1.,3.,5.,10.)
sens.par$frequency <- rev(c(2.,5.))#c(1.,2.,5.)
sens.par$t.res     <- rev(c(0.1,1.,2.,5.,10.)/365.25) 

# sensitivity run
k <- 0
while(k == 0){
for (i in 1:length(sens.par$frequency)){
  for (j in 1:length(sens.par$t.res)){
  
    if (((i == 2)&(j==2))) k <- 1
      PL$frequency <- sens.par$frequency[i]
      PL$t.res     <- sens.par$t.res[j]
      
      print(paste("start run: freq =",PL$frequency,"yr-1 t.res =",PL$t.res*365.25,"d"))
      
      PL <- set.parameter.list(PL)
      
      #  Experiment name
      Run.name <- paste("DredgeExp_Coastal_freq",sens.par$frequency[i],"tres",sens.par$t.res[j]*365.25,sep="_")
      
      ModelResults <- TrawlImpactRecovery(PL,SV.ini=SV.ini)
      
      plotProfile(model.natural=ModelResults$undisturbed.state,
                  model.chronic=ModelResults$disturbed.state$disturbed.state,
                  timepoint.chronic=NULL,
                  PL=ModelResults$PL,
                  ylimz = c(20., 0),
                  colors=magma(5) )
      
      save(ModelResults, file = paste("../Output/",Run.name,"_results.Rdata",sep=""))
      dev.off()
      #}    
  
}}
}
