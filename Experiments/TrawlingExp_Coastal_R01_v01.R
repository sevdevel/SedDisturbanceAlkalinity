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

# load baseline file
load("../Output/Baseline_DepthTransect.Rdata")
PL     <- ModelResults[[1]]$PL
SV.ini <- ModelResults[[1]]$output$y

# Adapt disturbance parameters 
PL$yrs <- 60
# PD in mud= 3.2 +/- 1. 95%CI[1.4, 5.2]
# Average erosion Beam Trawl
PL$mud.content <- 50.
#PL$depth.ero <- sum((2.602*50/100+1.206e-3*c(19,1595,135,699) + 1.321e-2*50/100*c(19,1595,135,699))/(2.6*1e3)*1e2)  # Beam trawl small vessel
#PL$depth.ero <- sum((2.602*50/100+1.206e-3*c(1158,7.6,4020) + 1.321e-2*50/100*c(1158,7.6,4020))/(2.6*1e3)*1e2)      # otter trawl 'average Scottish' vessel
PL$depth.mix  <- 1. 
PL$res.effect <- c(1.,1.,1.) 
PL$longevity   <- 7.5
PL$TI.Benthos  <- FALSE
# Update parameter list 
PL <- set.parameter.list(PL)

#create sensitivity range
sens.par <- list()
#sens.par$frequency <- c(1.,3.,5.,10.)
sens.par$frequency <- c(1.,2.,5.,10.)
sens.par$depth.mix <- c(2.2,3.3,4.4)
sens.par$t.res     <- c(1.5, 6.8, 30.)/365.25 
sens.par$depth.ero <- c(3.0,3.2,3.4)
sens.par$name      <- c("25","50","75")
# sensitivity run
# For uncertainty: take 25% and 75% percentile for resuspension time from Khedri et al.       -> 25%= 1.5, 50%= 6.8, 75%= 30
#                  take 25% and 75% percentile for sediment resuspension from Khedri et al.   -> 25%= 1.5, 50%= 1.6, 75%= 1.7
#                  take 25% and 75% percentile for sediment homogenization from Khedri et al. -> 25%= 2.2, 50%= 3.3, 75%= 4.4 (mud) 
#                                                                                             -> 25%= 1.3, 50%= 2.3, 75%= 3.4 (sand) 


for (i in 1:length(sens.par$frequency)){
  for (j in c(1,2,3)){
  
    #if (!((i == 1)&(j==1))){
      PL$frequency <- sens.par$frequency[i]
      PL$depth.mix <- sens.par$depth.mix[j]
      PL$t.res     <- sens.par$t.res[j]
      PL$depth.ero <- sens.par$depth.ero[j]
      
      print(paste("start run: freq =",PL$frequency,"yr-1 percentile =",sens.par$name[j]))
      
      PL <- set.parameter.list(PL)
      
      #  Experiment name
      Run.name <- paste("TrawlingExp_Coastal_freq",sens.par$frequency[i],"perc",sens.par$name[j],sep="_")
      
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

