#####################################################################################################
# Auxiliary functions to assess impact of trawling on benthic DIC/Alk balance 
# Author: Sebastiaan van de Velde - Royal Belgian Institute of Natural Sciences
#                                   Université Libre de Bruxelles
# contact: svandevelde@naturalsciences.be
#####################################################################################################

library(AquaEnv)

#-----------------------------------------------------------------------------
# Function: init.diffusion.coefficient 
#-----------------------------------------------------------------------------

init.diffusion.coefficient <- function (species,S,TC,P,conv.fac,grid,tort.grid)
{
  Dmol <- conv.fac*diffcoeff(S=S,t=TC,P=P,species=species)[[species]]
  D.grid <- setup.prop.1D(value=Dmol,grid=grid)
  D.grid$mid <- D.grid$mid/tort.grid$mid
  D.grid$int <- D.grid$int/tort.grid$int
  return(list(D.grid=D.grid,Dmol=Dmol))
}

#-----------------------------------------------------------------------------
# Function: Integrate.conc 
#-----------------------------------------------------------------------------

Integrate.conc <- function(C,type="Solid",PL,N.end=length(C)){
  
  if (type=="Solid"){
    Int.C <- sum((PL$grid$dx*PL$svf.grid$mid*C)[1:N.end])
  }
  if (type=="Liquid"){
    Int.C <- sum((PL$grid$dx*PL$por.grid$mid*C)[1:N.end])
  }
  return(Int.C)
}

#-----------------------------------------------------------------------------
# Function: Integrate.reac 
#-----------------------------------------------------------------------------

Integrate.reac <- function(R,PL,N.end=length(R)){
  
  Int.R <- sum((PL$grid$dx*R)[1:N.end])

  return(Int.R)
}

#====================================================================================================
# Perturbation functions
#====================================================================================================

#-----------------------------------------------------------------------------
# Function: Perturb.sed and Perturb.sediment
#-----------------------------------------------------------------------------

Perturb.sed <- function(C,PL,type="Solid",type.per=c("erode","mix"),
                        depth.mix=0.0,C.ow=NULL,particle.sorting=F,
                        depth.ero=0.0,t.susp=0.0,k.reac=0.0){
  
  new.C     <- C
  InstaFlux <- 0.
  ImpactDepth <- 0.
  
  if ("erode"%in%type.per){
    
    ImpactDepth <- depth.ero
    
    if (any(PL$grid$x.mid<=ImpactDepth)){
      N.ero <- max(which(PL$grid$x.mid<=ImpactDepth))
      
      if (type=="Solid"){
        Int.C          <- Integrate.conc(C=new.C,type="Solid",PL=PL,N.end=N.ero)
        new.C[1:N.ero] <- Decay.model(C.0=C[1:N.ero],k=k.reac,t=t.susp)$C.t
        Int.C.new      <- Integrate.conc(C=new.C,type="Solid",PL=PL,N.end=N.ero)
        
        InstaFlux <- InstaFlux + (Int.C.new-Int.C) # umol cm-2
      }
      if (type=="Liquid"){
        # nothing happens, as all impact is done in the mixing part
      }
    }
  }
  
  if ("mix"%in%type.per){
    
    ImpactDepth <- depth.mix
    
    if (any(PL$grid$x.mid<=ImpactDepth)){
      N.mix <- max(which(PL$grid$x.mid<=ImpactDepth))
      
      if (type=="Solid"){
        Int.C          <- Integrate.conc(C=new.C,type="Solid",PL=PL,N.end=N.mix)
        if (!(particle.sorting)){
          new.C[1:N.mix]     <- Int.C/sum((PL$grid$dx*PL$svf.grid$mid)[1:N.mix])
        } else {
          new.C[1:N.mix] <- (Int.C/length(1:N.mix))/(PL$grid$dx*PL$svf.grid$mid)[1:N.mix]
          # Introduce 'grain sorting' -> bigger particles sink faster, so fine ones end on top
        }
        Int.C.new      <- Integrate.conc(C=new.C,type="Solid",PL=PL,N.end=N.mix)
      }
      if (type=="Liquid"){
        Int.C          <- Integrate.conc(C=new.C,type="Liquid",PL=PL,N.end=N.mix)
        #new.C[1:(N.mix-5)]         <- C.ow
        #new.C[(N.mix-5):(N.mix+5)] <- new.C[(N.mix-5)] + (new.C[(N.mix+5)] - new.C[(N.mix-5)])/10*((N.mix-5):(N.mix+5)-((N.mix-5)))
        new.C[1:N.mix]         <- C.ow
        
        Int.C.new      <- Integrate.conc(C=new.C,type="Liquid",PL=PL,N.end=N.mix)
      }
      InstaFlux <- InstaFlux + (Int.C.new-Int.C) # umol cm-2
    }
  }
 
  return(list(new.C=new.C,InstaFlux=InstaFlux))
}

Perturb.sediment <- function(IniconsArray,InstaFluxArray,PL,particle.sorting=F){
  
  PerturbIniconsArray <- IniconsArray
  
  # Solids
  Perturb.temp <- Perturb.sed(IniconsArray[1:PL$N.gr],type="Solid",type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,particle.sorting=particle.sorting,
                              depth.ero=PL$depth.ero,t.susp=PL$t.res,k.reac=PL$k.f*PL$res.effect[1])
  
  PerturbIniconsArray[1:PL$N.gr] <- Perturb.temp$new.C
  InstaFluxArray$CH2O.f  <- InstaFluxArray$CH2O.f + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(PL$N.gr+1):(2*PL$N.gr)],type="Solid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,particle.sorting=particle.sorting,
                              depth.ero=PL$depth.ero,t.susp=PL$t.res,k.reac=PL$k.s.base*PL$res.effect[2])
  
  PerturbIniconsArray[(PL$N.gr+1):(2*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$CH2O.s              <- InstaFluxArray$CH2O.s + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(2*PL$N.gr+1):(3*PL$N.gr)],  type="Solid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,particle.sorting=particle.sorting,
                              depth.ero=PL$depth.ero,t.susp=PL$t.res,k.reac=PL$k.v.base*PL$res.effect[3])
  
  PerturbIniconsArray[(2*PL$N.gr+1):(3*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$CH2O.v                 <- InstaFluxArray$CH2O.v+ Perturb.temp$InstaFlux
  
 
  Perturb.temp <- Perturb.sed(IniconsArray[(12*PL$N.gr+1):(13*PL$N.gr)],type="Solid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,particle.sorting=particle.sorting,
                              depth.ero=PL$depth.ero,t.susp=PL$t.res,k.reac=PL$k.FeS2)
  
  PerturbIniconsArray[(12*PL$N.gr+1):(13*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$FeS2                     <- InstaFluxArray$FeS2 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(14*PL$N.gr+1):(15*PL$N.gr)],type="Solid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,particle.sorting=particle.sorting,
                              depth.ero=PL$depth.ero,t.susp=0.,k.reac=0.)
  
  PerturbIniconsArray[(14*PL$N.gr+1):(15*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$CaCO3                     <- InstaFluxArray$CaCO3 + Perturb.temp$InstaFlux
  
  
  Perturb.temp <- Perturb.sed(IniconsArray[(16*PL$N.gr+1):(17*PL$N.gr)],type="Solid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,particle.sorting=particle.sorting,
                              depth.ero=PL$depth.ero,t.susp=0.,k.reac=0.)
  
  PerturbIniconsArray[(16*PL$N.gr+1):(17*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$CaCO3.arg                  <- InstaFluxArray$CaCO3.arg + Perturb.temp$InstaFlux
  
 
  # Liquids
  Perturb.temp <- Perturb.sed(IniconsArray[(3*PL$N.gr+1):(4*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$O2.ow)
  PerturbIniconsArray[(3*PL$N.gr+1):(4*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$O2                  <- InstaFluxArray$O2 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(4*PL$N.gr+1):(5*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$NO3.ow)
  PerturbIniconsArray[(4*PL$N.gr+1):(5*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$NO3              <- InstaFluxArray$NO3 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(5*PL$N.gr+1):(6*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$SO4.ow)
  PerturbIniconsArray[(5*PL$N.gr+1):(6*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$SO4                 <- InstaFluxArray$SO4 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(6*PL$N.gr+1):(7*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$CO2.ow)
  PerturbIniconsArray[(6*PL$N.gr+1):(7*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$CO2                    <- InstaFluxArray$CO2 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(7*PL$N.gr+1):(8*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$NH4.ow)
  PerturbIniconsArray[(7*PL$N.gr+1):(8*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$NH4                    <- InstaFluxArray$NH4 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(8*PL$N.gr+1):(9*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$H2PO4.ow)
  PerturbIniconsArray[(8*PL$N.gr+1):(9*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$H2PO4              <- InstaFluxArray$H2PO4 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(9*PL$N.gr+1):(10*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$H2S.ow)
  PerturbIniconsArray[(9*PL$N.gr+1):(10*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$H2S                   <- InstaFluxArray$H2S + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(10*PL$N.gr+1):(11*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$N2.ow)
  PerturbIniconsArray[(10*PL$N.gr+1):(11*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$N2                     <- InstaFluxArray$N2 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(11*PL$N.gr+1):(12*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$TA.ow)
  PerturbIniconsArray[(11*PL$N.gr+1):(12*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$TA                     <- InstaFluxArray$TA + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(13*PL$N.gr+1):(14*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$CH4.ow)
  PerturbIniconsArray[(13*PL$N.gr+1):(14*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$CH4     <- InstaFluxArray$CH4 + Perturb.temp$InstaFlux
  
  Perturb.temp <- Perturb.sed(IniconsArray[(15*PL$N.gr+1):(16*PL$N.gr)],type="Liquid", type.per=c("erode","mix"),
                              PL=PL,depth.mix=PL$depth.mix,C.ow=PL$Ca.ow)
  PerturbIniconsArray[(15*PL$N.gr+1):(16*PL$N.gr)] <- Perturb.temp$new.C
  InstaFluxArray$Ca     <- InstaFluxArray$Ca + Perturb.temp$InstaFlux
  
  return(list(PerturbIniconsArray=PerturbIniconsArray,InstaFluxArray=InstaFluxArray))
}


#-----------------------------------------------------------------------------
# Function: decay.function 
#-----------------------------------------------------------------------------

Decay.model <- function(C.0,k,t){
  
  # model  
  
  C.t <- C.0*exp(-k*t)
  
  return(list(C.t=C.t,C.lost=(C.0-C.t)))
  
}

#-----------------------------------------------------------------------------
# Function: Benthos.Depletion and Benthos.Recovery
#-----------------------------------------------------------------------------

# Logistic growth modifier.
# This function simulates  logistic growth of a population with a certain 
# growth rate r, subject to impacts (at times impacts) which cull the population
# with a depletion factor d. This transience is scaled relative to the carrying
# capacity K.

Benthos.Recovery <- function(times, impacts, d, r) {
  K <- 1
  out <- NULL
  out <- cbind(1, out)
  
  for (i in 2:(length(times) + 1)) {
    N <- out[i - 1]
    if ( min(abs(times[i-1]-impacts))<1e-6 ) {
      newout <- pmax(0.1, N * pmax(0, (1 - d)))
    } else {
      newout <- N + r * N * (1 - N / K)
    }
    out <- c(out, newout)
  }
  return(out[-1])
}

# Calculates depletion of species abundance (d) based on the penetration 
# depth and mud concent according to Sciberras et al., (2018).

Benthos.Depletion <- function(p.depth, mud = 10) {
  d      <- p.depth * 3 + mud * 0.3
  d      <- d / 100
  return(round(d, 2))
}

#====================================================================================================
# Plotting functions
#====================================================================================================

# 6. plotProfile

plotProfile <- function(model.natural,
                        model.chronic=NULL,
                        model.recovery.nobiot=NULL,
                        model.recovery.biot=NULL,
                        timepoint.chronic=NULL,
                        PL,
                        ylimz = c(20., 0),
                        colors=magma(4)
) {
  
  Csolid.to.wtperc     <- 106*12.*1e-6/PL$rho.sed*1e2 # umol OM cm-3 solid    -> g C 100 g-1 solid
  CaCO3solid.to.wtperc <- 12.*1e-6/PL$rho.sed*1e2     # umol CaCO3 cm-3 solid -> g C 100 g-1 solid
  
  
  CH2O.f.natural <- model.natural$y[1:PL$N.gr]
  CH2O.s.natural <- model.natural$y[(PL$N.gr+1):(2*PL$N.gr)]
  CH2O.v.natural <- model.natural$y[(2*PL$N.gr+1):(3*PL$N.gr)]
  O2.natural     <- model.natural$y[(3*PL$N.gr+1):(4*PL$N.gr)]
  NO3.natural    <- model.natural$y[(4*PL$N.gr+1):(5*PL$N.gr)]
  SO4.natural    <- model.natural$y[(5*PL$N.gr+1):(6*PL$N.gr)]
  CO2.natural    <- model.natural$y[(6*PL$N.gr+1):(7*PL$N.gr)]
  NH4.natural    <- model.natural$y[(7*PL$N.gr+1):(8*PL$N.gr)]
  H2PO4.natural  <- model.natural$y[(8*PL$N.gr+1):(9*PL$N.gr)]
  H2S.natural    <- model.natural$y[(9*PL$N.gr+1):(10*PL$N.gr)]
  N2.natural     <- model.natural$y[(10*PL$N.gr+1):(11*PL$N.gr)]
  TA.natural     <- model.natural$y[(11*PL$N.gr+1):(12*PL$N.gr)]
  FeS2.natural   <- model.natural$y[(12*PL$N.gr+1):(13*PL$N.gr)]
  CH4.natural    <- model.natural$y[(13*PL$N.gr+1):(14*PL$N.gr)]
  CaCO3.natural  <- model.natural$y[(14*PL$N.gr+1):(15*PL$N.gr)]
  CaCO3.arg.natural  <- model.natural$y[(16*PL$N.gr+1):(17*PL$N.gr)]
  
  pH.natural <- model.natural$pH
  
  CH2O.lim  <- c(0.,max((CH2O.f.natural+CH2O.s.natural+CH2O.v.natural)*Csolid.to.wtperc))
  O2.lim    <- c(0.,PL$O2.ow*1e3)
  NO3.lim    <- c(0.,PL$NO3.ow*1e3)
  SO4.lim   <- c(0.,PL$SO4.ow)
  DIC.lim   <- c(0.,max(CO2.natural))
  TA.lim    <- c(0.,max(TA.natural))
  FeS2.lim  <- c(0.,max(FeS2.natural))
  CaCO3.lim <- c(min(c(CaCO3.natural,CaCO3.arg.natural)),max(c(CaCO3.natural,CaCO3.arg.natural)))*CaCO3solid.to.wtperc
  pH.lim <- range(pH.natural)
  
  if (!is.null(model.chronic)){
    
    if (is.null(timepoint.chronic)) timepoint.chronic <- nrow(model.chronic)
    Chronic.time <- model.chronic[timepoint.chronic,1]
    
    CH2O.f.chronic <- model.chronic[timepoint.chronic, 2            :  (PL$N.gr+1)] 
    CH2O.s.chronic <- model.chronic[timepoint.chronic, (PL$N.gr+2)  : (2*PL$N.gr+1)] 
    CH2O.v.chronic <- model.chronic[timepoint.chronic, (2*PL$N.gr+2): (3*PL$N.gr+1)] 
    O2.chronic     <- model.chronic[timepoint.chronic, (3*PL$N.gr+2): (4*PL$N.gr+1)] 
    NO3.chronic    <- model.chronic[timepoint.chronic, (4*PL$N.gr+2): (5*PL$N.gr+1)] 
    SO4.chronic    <- model.chronic[timepoint.chronic, (5*PL$N.gr+2): (6*PL$N.gr+1)] 
    CO2.chronic    <- model.chronic[timepoint.chronic, (6*PL$N.gr+2): (7*PL$N.gr+1)] 
    NH4.chronic    <- model.chronic[timepoint.chronic, (7*PL$N.gr+2): (8*PL$N.gr+1)] 
    H2PO4.chronic  <- model.chronic[timepoint.chronic, (8*PL$N.gr+2): (9*PL$N.gr+1)] 
    H2S.chronic    <- model.chronic[timepoint.chronic, (9*PL$N.gr+2):(10*PL$N.gr+1)] 
    N2.chronic     <- model.chronic[timepoint.chronic,(10*PL$N.gr+2):(11*PL$N.gr+1)] 
    TA.chronic     <- model.chronic[timepoint.chronic,(11*PL$N.gr+2):(12*PL$N.gr+1)] 
    FeS2.chronic   <- model.chronic[timepoint.chronic,(12*PL$N.gr+2):(13*PL$N.gr+1)] 
    CH4.chronic    <- model.chronic[timepoint.chronic,(13*PL$N.gr+2):(14*PL$N.gr+1)] 
    CaCO3.chronic  <- model.chronic[timepoint.chronic,(14*PL$N.gr+2):(15*PL$N.gr+1)] 
    CaCO3.arg.chronic  <- model.chronic[timepoint.chronic,(16*PL$N.gr+2):(17*PL$N.gr+1)] 
    
    pH.chronic  <- model.chronic[timepoint.chronic,(17*PL$N.gr+2):(18*PL$N.gr+1)] 
    
    CH2O.lim  <- c(0.,max(c(CH2O.lim,(CH2O.f.chronic+CH2O.s.chronic+CH2O.v.chronic)*Csolid.to.wtperc)))
    DIC.lim   <- c(0.,max(c(DIC.lim,CO2.chronic)))
    TA.lim    <- c(0.,max(c(TA.lim,TA.chronic)))
    CaCO3.lim <- c(min(c(CaCO3.lim,CaCO3.chronic*CaCO3solid.to.wtperc,CaCO3.arg.chronic*CaCO3solid.to.wtperc)),
                   max(c(CaCO3.lim,CaCO3.chronic*CaCO3solid.to.wtperc,CaCO3.arg.chronic*CaCO3solid.to.wtperc)))
    pH.lim <- range(c(pH.natural,pH.chronic))
  }
  if (!is.null(model.recovery.nobiot)){
  }
  if (!is.null(model.recovery.biot)){
  }
  
  
  par(mfrow=c(3,3),oma=c(0,2,4,1))
  
  plot(x=(CH2O.f.natural+CH2O.s.natural+CH2O.v.natural)*Csolid.to.wtperc,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       xlab="",ylab="",type='l',lwd=2,col=colors[1],xlim=CH2O.lim )
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("POC (wt%)", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=(CH2O.f.chronic+CH2O.s.chronic+CH2O.v.chronic)*Csolid.to.wtperc,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  plot(x=O2.natural*1e3,y=PL$grid$x.mid,ylim=c(5.,0.),axes=F,
       xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=O2.lim)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("O2 (µM)", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=O2.chronic*1e3,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  if (!is.null(model.chronic)){
    mtext(paste("Chronic Trawling for",Chronic.time,"yrs",sep=" "), side = 3, line=4.5, font=2, cex=1.5)
  }
  
  #plot(x=NO3.natural*1e3,y=PL$grid$x.mid,ylim=c(5.,0.),axes=F,
  #     xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=NO3.lim)
  #axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  #mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  #axis(side = 3, lwd = 2, cex.axis= 1.3)
  #mtext("NO3 (µM)", side = 3, line = 2.5, font = 2, cex = 1.1)
  #if (!is.null(model.chronic)){
  #  lines(x=NO3.chronic*1e3,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  #}

  plot(x=pH.natural,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=pH.lim)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("pH", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=pH.chronic,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  plot(x=SO4.natural,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=SO4.lim)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("SO4 (mM)", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=SO4.chronic,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  #plot(x=H2S.natural*1e3,y=PL$grid$x.mid,ylim=ylimz,axes=F,
  #     xlab="",ylab="",type='l',lwd=2, col=colors[1])
  #axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  #mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  #axis(side = 3, lwd = 2, cex.axis= 1.3)
  #mtext("H2S (µM)", side = 3, line = 2.5, font = 2, cex = 1.1)
  #if (!is.null(model.chronic)){
  #  lines(x=H2S.chronic*1e3,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  #}
  
  plot(x=FeS2.natural/PL$rho.sed,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=FeS2.lim/PL$rho.sed)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("FeS2 (umol g-1)", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=FeS2.chronic/PL$rho.sed,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  plot(x=CO2.natural,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=DIC.lim)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("DIC (mM)", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=CO2.chronic,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  plot(x=TA.natural,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       xlab="",ylab="",type='l',lwd=2, col=colors[1],xlim=TA.lim)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("TA (mM)", side = 3, line = 2.5, font = 2, cex = 1.1)
  if (!is.null(model.chronic)){
    lines(x=TA.chronic,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  }
  
  plot(x=CaCO3.natural*CaCO3solid.to.wtperc,y=PL$grid$x.mid,ylim=ylimz,
       axes=F,xlab="",ylab="",type='l',lwd=2, col=colors[1], xlim=CaCO3.lim)
  axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  axis(side = 3, lwd = 2, cex.axis= 1.3)
  mtext("CaCO3 (wt%)", side = 3, line = 2.5, font = 2, cex = 1.1)
  lines(x=CaCO3.arg.natural*CaCO3solid.to.wtperc,y=PL$grid$x.mid,lty=3,lwd=2,col=colors[1])
  if (!is.null(model.chronic)){
    lines(x=CaCO3.chronic*CaCO3solid.to.wtperc,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
    lines(x=CaCO3.arg.chronic*CaCO3solid.to.wtperc,y=PL$grid$x.mid,lty=4,lwd=2,col=colors[2])
  }
  
  #OhmC.natural <- carb(flag=15,var1=TA.natural*1e-3,var2=CO2.natural*1e-3)$OmegaCalcite
  #plot(x=OhmC.natural,y=PL$grid$x.mid,ylim=ylimz,axes=F,
       #     xlab="",ylab="",type='l',lwd=2, col=colors[1])
  #axis(side = 2, lwd = 2, cex.axis= 1.3, las = 1)
  #mtext("Depth (cm)", side = 2, line = 2.5, font = 2, cex = 1.1)
  #axis(side = 3, lwd = 2, cex.axis= 1.3)
  #mtext("Ohm Calcite", side = 3, line = 2.5, font = 2, cex = 1.1)
  #if (!is.null(model.chronic)){
    #OhmC.chronic <- carb(flag=15,var1=TA.chronic*1e-3,var2=CO2.chronic*1e-3)$OmegaCalcite
    #  lines(x=OhmC.chronic,y=PL$grid$x.mid,lty=2,lwd=2,col=colors[2])
  #}
  
  plot.new()
  legend("left",bty='n',legend=c("Natural state","Chronic Trawling"),lty=c(1,2),lwd=2,col=colors[1:2])
}

#====================================================================================================
# Parameter function
#====================================================================================================

set.parameter.list <- function(PL=NULL)
{

  if (is.null(PL)) PL <- list()
  
  # Number of state variables and reactions
  PL$N.var  <- 17 
  PL$N.reac <- 15 
  # Transient simulation parameters
  if(is.null(PL$duration)) PL$duration <- 15 # Duration of the transient simulation [yr]

  # Model domain and grid definition
  if(is.null(PL$N.gr)) PL$N.gr <- 200          # total number of grid layers
  if(is.null(PL$L))    PL$L    <- 20          # total depth of sediment domain [cm]
  if(is.null(PL$dx.1)) PL$dx.1 <- PL$L/(PL$N.gr*10) # size of the first grid cell [cm]
  PL$grid <- setup.grid.1D(x.up = 0, L = PL$L, N = PL$N.gr, dx.1 = PL$dx.1)
  
  if(is.null(PL$C.lim)) PL$C.lim <- 1E-2
    
  # Environmental parameters for the calculation of diffusion coefficients
  if(is.null(PL$S))           PL$S           <- 35    # salinity
  if(is.null(PL$TC))          PL$TC          <- 10    # temperature [deg C]
  if(is.null(PL$water.depth)) PL$water.depth <- 50. # water depth
  
  PL$P <- (PL$water.depth/10+1.)*1.013 # pressure [bar]
  
  # Sediment model parameters
  if(is.null(PL$rho.sed)) PL$rho.sed <- 2.6  # density solid sediment [g cm-3]
  
  if(is.null(PL$por.0))   PL$por.0   <- 0.8 # porosity at the SWI 
  if(is.null(PL$por.inf)) PL$por.inf <- 0.8 # porosity at infinite depth
  if(is.null(PL$por.att)) PL$por.att <- 3.  # porosity attenuation coefficient [cm]
  PL$por.grid <- setup.prop.1D(func=p.exp, grid=PL$grid, y.0 = PL$por.0, y.inf = PL$por.inf, x.att=PL$por.att)
  PL$svf.grid <- setup.prop.1D(value = 0, grid = PL$grid)
  PL$svf.grid$mid <- (1-PL$por.grid$mid)
  PL$svf.grid$int <- (1-PL$por.grid$int)
    
  PL$tort.grid     <- setup.prop.1D(value = 0.0, grid = PL$grid)
  PL$tort.grid$mid <- 1 - 2*log(PL$por.grid$mid)
  PL$tort.grid$int <- 1 - 2*log(PL$por.grid$int)
    
  if(is.null(PL$v.sed)) PL$v.sed <- 0.2 # advection rate [cm yr-1]
  PL$v.grid <- setup.compaction.1D(v.inf = PL$v.sed, por.0=PL$por.0,por.inf=PL$por.inf, por.grid = PL$por.grid)$v
  PL$u.grid <- setup.compaction.1D(v.inf = PL$v.sed, por.0=PL$por.0,por.inf=PL$por.inf, por.grid = PL$por.grid)$u
  
  if(is.null(PL$Db.0))  PL$Db.0  <- 5.  # biodiffusitivity intensity [cm2 yr-1]
  if(is.null(PL$x.L))   PL$x.L   <- 5.  # biomixing depth [cm]
  if(is.null(PL$x.att)) PL$x.att <- 1.  # biomixing depth [cm]
  PL$Db.grid <- setup.prop.1D(func = p.sig, grid = PL$grid, y.0 = PL$Db.0, y.inf = 0., x.L = PL$x.L, x.att = PL$x.att)
    
  if(is.null(PL$irr.0))   PL$irr.0   <- 0.  # irrigation rate at SWI [yr-1]
  if(is.null(PL$irr.inf)) PL$irr.inf <- 0.  # irrigation rate at infinite depth [yr-1]
  if(is.null(PL$x.irr))   PL$x.irr   <- 3.5 # attenuation depth [cm]
  PL$irr.grid <- setup.prop.1D(func=p.exp,y.0=PL$irr.0,y.inf=PL$irr.inf,x.att=PL$x.irr,grid=PL$grid)
  
  if(is.null(PL$irr.O2))    PL$irr.O2     <- 1.
  if(is.null(PL$irr.NO3))   PL$irr.NO3    <- 1.
  if(is.null(PL$irr.SO4))   PL$irr.SO4    <- 1.
  if(is.null(PL$irr.CO2))   PL$irr.CO2    <- 1.
  if(is.null(PL$irr.NH4))   PL$irr.NH4    <- 0.
  if(is.null(PL$irr.H2PO4)) PL$irr.H2PO4  <- 1.
  if(is.null(PL$irr.CH4))   PL$irr.CH4    <- 0.
  if(is.null(PL$irr.H2S))   PL$irr.H2S    <- 0.
  if(is.null(PL$irr.N2))    PL$irr.N2     <- 1.
  if(is.null(PL$irr.TA))    PL$irr.TA     <- 1.
  if(is.null(PL$irr.Ca))    PL$irr.Ca     <- 1.
  
  # Degradation rate constants
  if(is.null(PL$k.f)) PL$k.f <- 10.    # fast degrading OC fraction [yr-1]
  if(is.null(PL$k.s)) PL$k.s.base <- 0.1    # slow degrading OC fraction [yr-1]
  if(is.null(PL$k.v)) PL$k.v.base <- 0.0001 # very slow degrading OC fraction [yr-1]

  if(is.null(PL$k.O2.acc)) PL$k.O2.acc <- 10. # acceleration factor aerobic degredation (-)
    
  if(is.null(PL$UL.SP)) PL$UL.SP <- 1.5 # upper limit of self-priming effect (1.5*k)
  if(is.null(PL$SP.SP)) PL$SP.SP <- 0.25 # steepness factor of self-priming effect
  ### calculated as k + (1.5 - 1)*k * (1 - exp(-SP.SP*[CH2O.f]))
  
  if(is.null(PL$model)) PL$model <- "m1" 
  
  # Fractions of OC
  if(is.null(PL$p.f)) PL$p.f <- 0.6  # fast degrading OC fraction
  if(is.null(PL$p.s)) PL$p.s <- 0.3 # slow degrading OC fraction
  if(is.null(PL$p.v)) PL$p.v <- 0.1 # very slow degrading OC fraction 
  
  #if((PL$p.f+PL$p.s+PL$p.c)!=1.0){
  #  print(paste("WARNING: sum of organic matter fractions is ",(PL$p.f+PL$p.s+PL$p.c)," which is arguably not 1 ...",sep=""))
  #}
  # Upper boundary conditions
  if(is.null(PL$Cforcing)) PL$Cforcing <- "constant"
  if(is.null(PL$Cflux))    PL$Cflux    <- 20.*1e3*1e-4*365.25 # carbon deposition flux [umol C cm-2 yr-1]
  if(is.null(PL$F.FeS2))   PL$F.FeS2   <- 0.                  # pyrite deposition flux [umol C cm-2 yr-1]
  if(is.null(PL$F.CaCO3))  PL$F.CaCO3  <- 0.                  # carbonate deposition flux [umol C cm-2 yr-1]
  if(is.null(PL$F.CaCO3.arg)) PL$F.CaCO3.arg <- 0.            # aragonite deposition flux [umol C cm-2 yr-1]
  
  if(is.null(PL$O2.ow))    PL$O2.ow    <- 0.35 # O2 concentration in overlaying water [umol cm-3]
  if(is.null(PL$NO3.ow))   PL$NO3.ow   <- 0.05 # NO3 concentration in overlaying water [umol cm-3]
  if(is.null(PL$NH4.ow))   PL$NH4.ow   <- 0.0  # NH4 concentration in overlaying water [umol cm-3]
  if(is.null(PL$H2PO4.ow)) PL$H2PO4.ow <- 0.0  # H2PO4 concentration in overlaying water [umol cm-3]
  if(is.null(PL$N2.ow))    PL$N2.ow    <- 0.0  # N2 concentration in overlaying water [umol cm-3]
  if(is.null(PL$CO2.ow))   PL$CO2.ow   <- 2.1  # DIC concentration in overlaying water [umol cm-3]
  if(is.null(PL$TA.ow))    PL$TA.ow    <- 2.2  # AT concentration in overlaying water [umol cm-3]
  if(is.null(PL$SO4.ow))   PL$SO4.ow   <- 28.  # SO4 concentration in overlaying water [umol cm-3]
  if(is.null(PL$H2S.ow))   PL$H2S.ow   <- 0.   # H2S concentration in overlaying water [umol cm-3]
  if(is.null(PL$CH4.ow))   PL$CH4.ow   <- 0.   # CH4 concentration in overlaying water [umol cm-3]
  if(is.null(PL$Ca.ow))    PL$Ca.ow    <- 10.  # Ca concentration in overlaying water [umol cm-3]
  
  if(is.null(PL$O2.ds))    PL$O2.ds    <- NA     # concentration O2 at bottom boundary  [?mol cm-3]
  if(is.null(PL$NO3.ds))   PL$NO3.ds   <- NA     # concentration NO3 at bottom boundary  [?mol cm-3]
  if(is.null(PL$SO4.ds))   PL$SO4.ds   <- NA     # concentration SO4 at bottom boundary [?mol cm-3]
  if(is.null(PL$CO2.ds))   PL$CO2.ds   <- NA     # concentration DIC at bottom boundary [?mol cm-3]
  if(is.null(PL$NH4.ds))   PL$NH4.ds   <- NA     # concentration NH4 at bottom boundary [?mol cm-3]
  if(is.null(PL$H2PO4.ds)) PL$H2PO4.ds <- NA     # concentration H2PO4 at bottom boundary  [?mol cm-3]
  if(is.null(PL$H2S.ds))   PL$H2S.ds   <- NA     # concentration H2S at bottom boundary  [?mol cm-3]
  if(is.null(PL$N2.ds))    PL$N2.ds    <- NA     # concentration N2 at bottom boundary  [?mol cm-3]
  if(is.null(PL$TA.ds))    PL$TA.ds    <- NA     # concentration AT at bottom boundary  [?mol cm-3]
  if(is.null(PL$CH4.ds))   PL$CH4.ds   <- NA     # concentration CH4 at bottom boundary  [?mol cm-3]
  if(is.null(PL$Ca.ds))    PL$Ca.ds    <- NA     # concentration Ca at bottom boundary  [?mol cm-3]
  
  # Porosity and solid volume fraction profiles
    
  # Transport parameters

  # Diffusion coefficients calculated using the marelac package [m2 s-1]
  conv.fac <- 1e4*(3600 * 24 * 365.25) # [m2 s-1] to [cm2 yr-1]
  PL$D.O2.grid  <- init.diffusion.coefficient(species="O2",    S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.NO3.grid   <- init.diffusion.coefficient(species="NO3", S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.NH4.grid   <- init.diffusion.coefficient(species="NH4", S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.H2PO4.grid <- init.diffusion.coefficient(species="PO4", S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.N2.grid    <- init.diffusion.coefficient(species="N2",  S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.CO2.grid   <- init.diffusion.coefficient(species="HCO3",S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.TA.grid    <- init.diffusion.coefficient(species="HCO3",S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.SO4.grid   <- init.diffusion.coefficient(species="SO4", S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.H2S.grid   <- init.diffusion.coefficient(species="H2S", S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.CH4.grid   <- init.diffusion.coefficient(species="CH4", S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  PL$D.Ca.grid    <- init.diffusion.coefficient(species="Ca",  S=PL$S,TC=PL$TC,P=PL$P,conv.fac=conv.fac,grid=PL$grid,tort.grid=PL$tort.grid)$D.grid
  
  # Reaction parameters
  
  # OM mineralization reactions
  if(is.null(PL$Ks.O2))  PL$Ks.O2  <- 0.005 # [umol cm-3] Monod constant O2 reduction (Aerobic respiration)
  if(is.null(PL$Ks.NO3)) PL$Ks.NO3 <- 0.01  # [umol cm-3] Monod constant NO3 reduction (NO3 reduction)
  if(is.null(PL$Ks.SO4)) PL$Ks.SO4 <- 0.5   # [umol cm-3] Monod constant SO4 reduction (SO4 reduction)
  
  # Secondary redox reactions
  if(is.null(PL$k.Nit))  PL$k.Nit  <- 1e4 # [umol cm-3 yr-1] ammonium oxidation constant
  if(is.null(PL$k.CSO))  PL$k.CSO  <- 1e6 # [umol cm-3 yr-1] Canonical Sulfide Oxidation constant
  if(is.null(PL$k.AeOM)) PL$k.AeOM <- 1e4 # [umol cm-3 yr-1] rate constant of Aerobic Oxidation of Methane
  if(is.null(PL$k.AOM))  PL$k.AOM  <- 1e2 # [umol cm-3 yr-1] rate constant of Anaerobic Oxidation of Methane by sulfate
  
  if(is.null(PL$p_fe))   PL$p_fe   <- 0.1 # [-] fraction of sulphide that is buried as pyrite
  if(is.null(PL$k.FeS2)) PL$k.FeS2 <- 1e2 # [umol cm-3 yr-1] Pyrite oxidation rate constant
  if(is.null(PL$k.CaP))  PL$k.CaP  <- 0.4            # [umol cm-3 yr-1] calcite precipitation rate constant (Sulpis et al., 2022)
  if(is.null(PL$n.CaP))  PL$n.CaP  <- 1.76           # [-] calcite precipitation rate exponent (Zuddas and Mucci, 1998)
  if(is.null(PL$k.CaD))  PL$k.CaD  <- c(20.,6.3e-3)  # [yr-1] calcite dissolution rate constant (Sulpis et al., 2022)
  if(is.null(PL$n.CaD))  PL$n.CaD  <- c(4.7,0.11)    # [-] calcite dissolution rate exponent (Sulpis et al., 2022)
  if(is.null(PL$k.ArD))  PL$k.ArD  <- c(4.2e-2,1e-3)*1e2 # [yr-1] aragonite dissolution rate constant (Sulpis et al., 2022) -> aragonite never got calibrated, only calcite, so just revert to the originial one
  if(is.null(PL$n.ArD))  PL$n.ArD  <- c(1.46,0.13)   # [-] aragonite dissolution rate exponent (Sulpis et al., 2022)
  if(is.null(PL$OhmCa.crit))PL$OhmCa.crit  <- 0.828  # [-] critical calcite saturation (Sulpis et al., 2022)
  if(is.null(PL$OhmAr.crit))PL$OhmAr.crit  <- 0.835  # [-] critical aragonite saturation (Sulpis et al., 2022)
  
  #if(is.null(PL$Ksp.Cal))PL$Ksp.Cal  <- aquaenv(S=PL$S,t=PL$TC,Pa=PL$P)$Ksp_calcite   # [M^2] Calcite solubility product
  #if(is.null(PL$Ksp.Ara))PL$Ksp.Ara  <- aquaenv(S=PL$S,t=PL$TC,Pa=PL$P)$Ksp_aragonite # [M^2] Aragonite solubility product
  PL$Ksp.Cal  <- aquaenv(S=PL$S,t=PL$TC,Pa=PL$P)$Ksp_calcite   # [M^2] Calcite solubility product
  PL$Ksp.Ara  <- aquaenv(S=PL$S,t=PL$TC,Pa=PL$P)$Ksp_aragonite # [M^2] Aragonite solubility product
  
  # Perturbation parameters
  if(is.null(PL$frequency)) PL$frequency  <- 5          # [yr-1] frequency of disturbances 
  if(is.null(PL$yrs))       PL$yrs        <- 100        # [yr] how long trawling has been going on
  if(is.null(PL$t.res))     PL$t.res      <- 0.1/365.25 # [yr] time of resuspended particles to remain in solution
  if(is.null(PL$res.effect))PL$res.effect <- c(1.,1.,1.)# [-] multiplying factor for C mineralisation in suspension (f, s, v) fractions
  
  if(is.null(PL$TI.Benthos)) PL$TI.Benthos <- FALSE      # Logical - does trawling affect benthos?
  if(is.null(PL$Db.modif))   PL$Db.modif   <- cbind(c(1.,100.),c(1.,1.)) # [-] multiplying factor for bioturbation
  if(is.null(PL$mud.content))PL$mud.content<- 10.        # [%] mud content to calculated benthos depletion
  if(is.null(PL$longevity))  PL$longevity  <- 1.         # [yr] longevity of benthos
  
  if(is.null(PL$depth.mix)) PL$depth.mix  <- 2.2        # [cm] depth of mixed layer
  if(is.null(PL$DC.trawl))  PL$DC.trawl   <- 5e3         # [N m-1] Hydrodynamic drag of dragged gear
  
  if(is.null(PL$depth.ero)) {
    PL$depth.ero  <- (2.602*PL$mud.content/100+1.206e-3*PL$DC.trawl + 1.321e-2*PL$mud.content/100*PL$DC.trawl)/(PL$rho.sed*1e3)*1e2
  }  # [cm] depth of eroded layer, based on density and O'Neill and Ivanovic - link between resuspended sediment, mud content and hydrodynamic drag
  if(is.null(PL$particle.sorting)) PL$particle.sorting <- FALSE # logical - does trawling lead to particle sorting? 
  # Assembling the parameters list
  
  return (PL)
}

#====================================================================================================
# Set initial state
#====================================================================================================

set.initial.state <- function(PL,type="crude",previous.run=NULL){
  
  if (type=="crude"){
    
    # Set crude initial conditions
    
    CH2O.f.in <- rep(0,length.out = PL$N.gr)
    CH2O.s.in <- rep(0,length.out = PL$N.gr)
    CH2O.v.in <- rep(0,length.out = PL$N.gr)
    O2.in    <- rep(PL$O2.ow,length.out = PL$N.gr)
    NO3.in   <- rep(PL$NO3.ow,length.out = PL$N.gr)
    SO4.in   <- rep(PL$SO4.ow,length.out = PL$N.gr)
    CO2.in   <- rep(PL$CO2.ow,length.out = PL$N.gr)
    NH4.in   <- rep(PL$NH4.ow,length.out = PL$N.gr)
    H2PO4.in <- rep(PL$H2PO4.ow,length.out = PL$N.gr)
    H2S.in   <- rep(PL$H2S.ow,length.out = PL$N.gr)
    N2.in    <- rep(PL$N2.ow,length.out = PL$N.gr)
    TA.in    <- rep(PL$TA.ow,length.out = PL$N.gr)
    FeS2.in  <- rep(0.,length.out = PL$N.gr)
    CH4.in   <- rep(PL$CH4.ow,length.out = PL$N.gr)
    CaCO3.in <- rep(0.,length.out = PL$N.gr)
    Ca.in    <- rep(PL$Ca.ow,length.out = PL$N.gr)
    CaCO3.arg.in <- rep(0.,length.out = PL$N.gr)
    
  }
  
  if (type=="transient"){
    
    CH2O.f.in <- previous.run[nrow(previous.run), 2            :  (PL$N.gr+1)] 
    CH2O.s.in <- previous.run[nrow(previous.run), (PL$N.gr+2)  : (2*PL$N.gr+1)] 
    CH2O.v.in <- previous.run[nrow(previous.run), (2*PL$N.gr+2): (3*PL$N.gr+1)] 
    O2.in     <- previous.run[nrow(previous.run), (3*PL$N.gr+2): (4*PL$N.gr+1)] 
    NO3.in    <- previous.run[nrow(previous.run), (4*PL$N.gr+2): (5*PL$N.gr+1)] 
    SO4.in    <- previous.run[nrow(previous.run), (5*PL$N.gr+2): (6*PL$N.gr+1)] 
    CO2.in    <- previous.run[nrow(previous.run), (6*PL$N.gr+2): (7*PL$N.gr+1)] 
    NH4.in    <- previous.run[nrow(previous.run), (7*PL$N.gr+2): (8*PL$N.gr+1)] 
    H2PO4.in  <- previous.run[nrow(previous.run), (8*PL$N.gr+2): (9*PL$N.gr+1)] 
    H2S.in    <- previous.run[nrow(previous.run), (9*PL$N.gr+2):(10*PL$N.gr+1)] 
    N2.in     <- previous.run[nrow(previous.run),(10*PL$N.gr+2):(11*PL$N.gr+1)] 
    TA.in     <- previous.run[nrow(previous.run),(11*PL$N.gr+2):(12*PL$N.gr+1)] 
    FeS2.in   <- previous.run[nrow(previous.run),(12*PL$N.gr+2):(13*PL$N.gr+1)] 
    CH4.in    <- previous.run[nrow(previous.run),(13*PL$N.gr+2):(14*PL$N.gr+1)] 
    CaCO3.in  <- previous.run[nrow(previous.run),(14*PL$N.gr+2):(15*PL$N.gr+1)] 
    Ca.in     <- previous.run[nrow(previous.run),(15*PL$N.gr+2):(16*PL$N.gr+1)] 
    CaCO3.arg.in  <- previous.run[nrow(previous.run),(16*PL$N.gr+2):(17*PL$N.gr+1)] 
    
    }
    
    # Initialization state variables vector 
    
    SV.ini <- c(CH2O.f.in,CH2O.s.in,CH2O.v.in, O2.in, NO3.in, SO4.in, CO2.in, NH4.in, H2PO4.in, H2S.in, N2.in, TA.in, FeS2.in, CH4.in, CaCO3.in, Ca.in, CaCO3.arg.in)
    
    return(SV.ini)
}

#====================================================================================================
# Set auxiliary Array
#====================================================================================================

Set.auxiliary.arrays <- function(PL){
  
  Stoichiometry.array <- as.data.frame(matrix(rep(0.0,PL$N.var*PL$N.reac),nrow=PL$N.reac,ncol=PL$N.var))
  colnames(Stoichiometry.array) <- c("CH2O.f","CH2O.s","CH2O.v","O2","NO3","SO4","CO2","NH4","H2PO4",
                                     "CH4","H2S","N2","TA","FeS2","CaCO3","Ca")
  rownames(Stoichiometry.array) <- c("Cmin.f","Cmin.s","Cmin.v","AR","DN","SR","MG","Nit","CSO","AeOM","PyO","AOM","CaD","CaP")
 #                                    CH2O.f, CH2O.s, CH2O.v, O2,     NO3,                                  SO4, CO2,  NH4, H2PO4, CH4,                                  H2S,       N2,                                 TA,                           FeS2, CaCO3, Ca
  Stoichiometry.array["Cmin.f",] <- c(-1.,    0.,     0.,     0.,       0.,                                  0.,   0.,   0.,   0.,    0.,                                 0.,      0.,                                  0.,                            0.,      0.,  0.  )
  Stoichiometry.array["Cmin.s",] <- c(0.,    -1.,     0.,     0.,       0.,                                  0.,   0.,   0.,   0.,    0.,                                 0.,      0.,                                  0.,                            0.,      0.,  0.  )
  Stoichiometry.array["Cmin.v",] <- c(0.,     0.,    -1.,     0.,       0.,                                  0.,   0.,   0.,   0.,    0.,                                 0.,      0.,                                  0.,                            0.,      0.,  0.  )
  Stoichiometry.array["AR",]     <- c(0.,     0.,     0.,  -106.,       0.,                                  0., 106.,  16.,   1.,    0.,                                 0.,      0.,                                 15.,                            0.,      0.,  0.  )
  Stoichiometry.array["DN",]     <- c(0.,     0.,     0.,     0., -424./5.,                                  0., 106.,  16.,   1.,    0.,                                 0., 212./5.,                             499./5.,                            0.,      0.,  0.  )
  Stoichiometry.array["SR",]     <- c(0.,     0.,     0.,     0.,       0.,  -(848-742*PL$p_fe)/(16-15*PL$p_fe), 106.,  16.,   1.,    0., (848-1590*PL$p_fe)/(16-15*PL$p_fe),      0., (1936-1709*PL$p_fe)/(16-15*PL$p_fe), (424*PL$p_fe)/(16-15*PL$p_fe),      0.,  0.  )
  Stoichiometry.array["MG",]     <- c(0.,     0.,     0.,     0.,       0.,                                  0.,  58.,  16.,   1.,   58.,                                 0.,      0.,                                 15.,                            0.,      0.,  0.  )
  Stoichiometry.array["Nit",]    <- c(0.,     0.,     0.,    -2.,       1.,                                  0.,  58.,  -1.,   0.,    0.,                                 0.,      0.,                                 -2.,                            0.,      0.,  0.  )
  Stoichiometry.array["CSO",]    <- c(0.,     0.,     0.,    -2.,       0.,                                  1.,  58.,   0.,   0.,    0.,                                -1.,      0.,                                 -2.,                            0.,      0.,  0.  )
  Stoichiometry.array["AeOM",]   <- c(0.,     0.,     0.,    -2.,       0.,                                  0.,   1.,   0.,   0.,   -1.,                                 0.,      0.,                                  0.,                            0.,      0.,  0.  )
  Stoichiometry.array["PyO",]    <- c(0.,     0.,     0.,-15./4.,       0.,                                  2.,   0.,   0.,   0.,   -1.,                                 0.,      0.,                                 -4.,                           -1.,      0.,  0.  )
  Stoichiometry.array["AOM",]    <- c(0.,     0.,     0.,     0.,       0.,                                 -1.,   1.,   0.,   0.,   -1.,                                 1.,      0.,                                  2.,                            0.,      0.,  0.  )
  Stoichiometry.array["CaD",]    <- c(0.,     0.,     0.,     0.,       0.,                                  0.,   1.,   0.,   0.,    0.,                                 0.,      0.,                                  2.,                            0.,     -1.,  1.  )
  Stoichiometry.array["CaP",]    <- c(0.,     0.,     0.,     0.,       0.,                                  0.,  -1.,   0.,   0.,    0.,                                 0.,      0.,                                 -2.,                            0.,      1., -1.  )
  
  return(Stoichiometry.array)
}

#====================================================================================================
