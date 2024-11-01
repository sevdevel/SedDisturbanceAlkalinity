#####################################################################################################
# Script to simulate impact of trawling (mixing/erosion) on a sediment column
# Author: Sebastiaan van de Velde - Royal Belgian Institute of Natural Sciences
#                                   Université Libre de Bruxelles
# contact: svandevelde@naturalsciences.be
#####################################################################################################

## Load the package for use
library(ReacTran)
library(marelac)
library(viridis)

TrawlImpactRecovery <- function(PL,SV.ini){
  
  InstaFlux <- data.frame(CH2O.f=rep(0.0,PL$frequency*PL$yrs),CH2O.s=rep(0.0,PL$frequency*PL$yrs),CH2O.v=rep(0.0,PL$frequency*PL$yrs),
                          FeS2=rep(0.0,PL$frequency*PL$yrs),O2=rep(0.0,PL$frequency*PL$yrs),NO3=rep(0.0,PL$frequency*PL$yrs),
                          SO4=rep(0.0,PL$frequency*PL$yrs),CO2=rep(0.0,PL$frequency*PL$yrs),NH4=rep(0.0,PL$frequency*PL$yrs),
                          H2PO4=rep(0.0,PL$frequency*PL$yrs),H2S=rep(0.0,PL$frequency*PL$yrs),N2=rep(0.0,PL$frequency*PL$yrs),
                          TA=rep(0.0,PL$frequency*PL$yrs),CH4=rep(0.0,PL$frequency*PL$yrs),CaCO3=rep(0.0,PL$frequency*PL$yrs),
                          Ca=rep(0.0,PL$frequency*PL$yrs),CaCO3.arg=rep(0.0,PL$frequency*PL$yrs))
  
  #PL <- set.parameter.list(PL=PL)
  
  #SV.ini <- set.initial.state(PL)
  
  undisturbed.state <- steady.1D(y=SV.ini, func=DiagC.model, parms=PL, nspec=PL$N.var, pos=TRUE)
  print("Predisturbed state calculated")
  x11()
  plotProfile(model.natural=undisturbed.state,PL=PL)
  
  #print(undisturbed.state$F.CaCO3.up)
  #print(undisturbed.state$F.CaCO3.down)
  #print(undisturbed.state$F.CaCO3.arg.up)
  #print(undisturbed.state$F.CaCO3.arg.down)
  
  inicons <- undisturbed.state$y
  
  # Equally distribute disturbances throughout the year
  if (PL$frequency > 1){
    ts <- 1. + floor(365/(PL$frequency+1))*1:PL$frequency
  } else {
    if (PL$frequency == 1){
      ts <- 364./2
    }  else { ts <- 0.}
  }
  ts <- c(ts, ts + rep(seq(from = 365, to = 365 * (PL$yrs-1), by = 365), each = length(ts)))
  if (ts[1] != 0.) ts <- c(0.,round(ts/365,2))
  
  # Apply disturbance on benthos
  
  if (PL$TI.Benthos==TRUE){
    Depl.Factor <- Benthos.Depletion(p.depth=PL$depth.mix, mud=PL$mud.content)
    PL$Db.modif <- cbind(seq(0.,PL$yrs,by=0.01),
                         Benthos.Recovery(times=seq(0.,PL$yrs,by=0.01), impacts=ts[-1], d=Depl.Factor, r=1/PL$longevity))
    par(mfrow=c(1,1))
    plot(x=PL$Db.modif[,1],y=PL$Db.modif[,2],xlim=c(0,3),xlab="time (yr)",ylab="Db modifier",type="l")
    
  } else {
    PL$Db.modif <- 1.
  }
  
  
  # Run transient with disturbance events
  
  for (i in 2:length(ts)){
    
    # Run until disturbance
    
    times <- seq(ts[i-1],ts[i],by=0.01)
    #print(times)
    
    if (i==2){
      disturbed.state <- ode.1D(y=inicons, times=times, func=DiagC.model, parms=PL, nspec=PL$N.var, 
                                method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-10, rtol=1E-7, verbose=F)
      
    } else {
      transient.dist <- ode.1D(y=inicons, times=times, func=DiagC.model, parms=PL, nspec=PL$N.var, 
                               method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-10, rtol=1E-7, verbose=F)
      disturbed.state <- rbind(disturbed.state,transient.dist)
    }
    
    plotProfile(model.natural=undisturbed.state,
                model.chronic=disturbed.state,
                PL=PL)
    
    # Set new initial state
    
    new.ini <- set.initial.state(PL=PL,type="transient",previous.run=disturbed.state)
    
    # Apply disturbance on sediment
    
    Perturbation.list <- Perturb.sediment(IniconsArray=new.ini,InstaFluxArray=InstaFlux[(i-1),],PL=PL,particle.sorting=PL$particle.sorting)
    inicons           <- Perturbation.list$PerturbIniconsArray
    InstaFlux[(i-1),] <- Perturbation.list$InstaFluxArray
    
    print("perturbed initial state")
    
    print(paste((i-1),"th disturbance event of ",length(ts)-1," done",sep=""))
  }
  
  # When last disturbance has passed, run until final steady state
  
  times <- seq(ts[length(ts)],PL$yrs,by=0.01)
  
  transient.dist <- ode.1D(y=inicons, times=times, func=DiagC.model, parms=PL, nspec=PL$N.var, 
                           method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-10, rtol=1E-7, verbose=F)
  disturbed.state <- rbind(disturbed.state,transient.dist)
  
  print("Finished simulations")
  
  # return results
  
  return(list("undisturbed.state"=undisturbed.state,
              "disturbed.state"=list("disturbed.state"=disturbed.state,"Instant.Fluxes"=InstaFlux),
              "PL"=PL))
}
