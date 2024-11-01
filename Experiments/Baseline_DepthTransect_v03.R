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

source('../Functions/DiagCModel_v04.R')
source('../Functions/SolveSaphePackage.R')
source('../Functions/AuxFunctions_v04.R')
load("../../../FluxDatabase/ATFluxDatabase.Rdata")
load("../../../FluxDatabase/Fluxsummary.Rdata")

# Create parameter list 
PL <- set.parameter.list()

#==================================================
# Adapt boundary conditions per water depth
#==================================================

var <- "TA"
cols <- viridis(10)

#--------------------------------------------------
# Mud
#--------------------------------------------------

# data figure

par(mfrow=c(1,1))
plot(x=10.,y=Depth.Matrix["F_TA","Coastal"],log='x',pch=16,col='black',xlim=c(1,10000),ylim=c(0,20.),
     xlab="Water depth (m)",ylab="AT flux (meq m-2 d-1)",cex.axis=1.3,cex.lab=1.3)
legend("topright",bty='n',legend=c("data","model"),pch=c(16,15),col=c("black","red"),cex=1.3)
points(x=100.,y=Depth.Matrix["F_TA","Shelf"],pch=16,col="black")
points(x=1000.,y=Depth.Matrix["F_TA","Slope"],pch=16,col="black")
points(x=5000.,y=Depth.Matrix["F_TA","Deep-sea"],pch=16,col="black")

# model runs

ModelResults <- list()
water.depth <- c(10.,100.,1000.,5000.) # m
TA.BW  <- c(2.3,2.3,2.34,2.38)   # from ocean biogeochemical dynamics textbook
DIC.BW <- c(2.02,2.12,2.26,2.28) # from ocean biogeochemical dynamics textbook

for (i in 1:4) print(pH.equilibration(S=35,T=c(15.15,10.29,2,2)[i],P=(1+water.depth[i]/10.)*1.013,TA=TA.BW[i]*1e-3,SumCO2=DIC.BW[i]*1e-3)$pH)
#for (i in 1:length(water.depth)){

PL <- set.parameter.list(PL=NULL)

PL$k.f <- 1.*10.  # yr-1 10x higher because i don't think 1 is high enough
PL$k.s <- 0.1 # yr-1
PL$p.f <- 0.7
PL$p.s <- 0.15
PL$p.v <- 0.15
PL$k.Nit  <- 1.e4 # umol cm-3 yr-1
PL$k.CSO  <- 1.e6 # umol cm-3 yr-1
PL$k.FeS2 <- 10.  # umol cm-3 yr-1
PL$x.irr <- 2. # Dale et al., 2015 - shelf irrigation
PL$p_fe        <- 0.5

PL$n.CaD <- c(4.7,0.11)       # Based on Naviaux et al. for 12°C
PL$k.CaD <- c(20.,6.3e-3)    #*30 Based on Naviaux et al. for 12°C (just increased the Sulpis value to a higher one)
PL$k.ArD <-  c(4.2e-2,3.8e-3)*1e2
PL$n.ArD  <- c(1.46,0.13)
PL$k.CaP <- 0.4 

PL$irr.H2S <- 0.6
PL$irr.CO2 <- 0.2
PL$irr.NH4 <- 0.3

  i <- 1
  
  if (i == 1){
    PL$L    <- 200
    PL$N.gr <- 2000
  } else{
    PL$L    <- 20
    PL$N.gr <- 200
  }
  
  PL$water.depth <- water.depth[i] 
  PL$irr.0       <- 0.8*365.25*(5.2 * 10^(0.76241122-0.00039724*PL$water.depth))/30.# assuming irrigation of 1 d-1 is highest possible and irrigation scales with bioturbation
  PL$Db.0        <- 5.2 * 10^(0.76241122-0.00039724*PL$water.depth) 
  PL$x.L         <- 1.+9.*(1.-exp(-PL$Db.0/3.))
  PL$v.sed       <- 3.3 * 10^(-0.87478367-0.00043512*PL$water.depth) 
  PL$Cflux       <- 0.53*40.61*PL$water.depth^(-0.571)*1e6*1e-4 # umol cm-2 yr-1
  #PL$F.CaCO3     <- 2./3.*0.06*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
  #PL$F.CaCO3.arg <- 1./3.*0.06*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
  PL$F.CaCO3     <- 2./3.*0.3*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
  PL$F.CaCO3.arg <- 1./3.*0.3*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
  PL$NO3.ow      <- 0.003734 + 0.0000584*PL$water.depth
  PL$CO2.ow      <- DIC.BW[i]
  PL$TA.ow       <- TA.BW[i]
  PL$TC          <- 15.69 - 0.054*PL$water.depth
  if(PL$TC<2.) PL$TC <- 2.
  
  # Update parameter list 
  PL <- set.parameter.list(PL)
  SV.ini <- set.initial.state(PL)
  #if (i == 1) SV.ini <- set.initial.state(PL)
  #if (i != 1) SV.ini <- ModelResults[[i-1]]$output$y
    
  PL$Run.name <- paste("Baseline_DepthTransect_withirr",PL$water.depth,sep="_")
  
  out.temp <- steady.1D(y=SV.ini, func=DiagC.model, parms=PL, nspec=PL$N.var, pos=TRUE)
  print(paste("water depth=",PL$water.depth,"m"))
  print(paste("AT flux=",-((out.temp$F.TA.up+out.temp$TA.irr.int) - (out.temp$F.H2S.up+out.temp$H2S.irr.int))*1e-3*1e4/365.25,"meq m-2 d-1"))
  print(paste("DIC flux=",-((out.temp$F.CO2.up+out.temp$CO2.irr.int))*1e-3*1e4/365.25,"mmol m-2 d-1"))
  print(paste("O2 flux=",((out.temp$F.O2.up+out.temp$O2.irr.int) - (out.temp$F.H2S.up+out.temp$H2S.irr.int))*1e-3*1e4/365.25,"mmol m-2 d-1"))
  
  if (i == 1) print(Depth.Matrix[,"Coastal"])
  if (i == 2) print(Depth.Matrix[,"Shelf"])
  if (i == 3) print(Depth.Matrix[,"Slope"])
  if (i == 4) print(Depth.Matrix[,"Deep-sea"])
  
  x11()
  plotProfile(model.natural=out.temp,
              model.chronic=NULL,
              PL=PL,
              ylimz = c(200., 0))
  
  SS <- attributes(out.temp)$steady
  if(!SS){
    SV.ini <- out.temp$y
    times <- c(0.,0.5,1.,5.,10.,100.,500.,1000.,5000.)
    transient <- ode.1D(times=times, y=SV.ini, func=DiagC.model, parms=PL, nspec=PL$N.var, method="vode",maxsteps=100000,verbose=TRUE)
    #if (i == 4) {
    #  SV.ini <- transient[nrow(transient),1:(PL$N.gr*PL$N.var)]
    #  if (i == 4) times <- c(0.,0.5,1.,3.)
    #  transient <- ode.1D(times=times, y=SV.ini, func=DiagC.model, parms=PL, nspec=PL$N.var, method="vode",maxsteps=100000,verbose=TRUE)
    #}
    ModelResults[[i]] <- list(PL=PL,output=transient[nrow(transient),])
  } else{
    ModelResults[[i]] <- list(PL=PL,output=out.temp)
    }

  
save(ModelResults,file="Baseline_DepthTransect.Rdata")
  
par(mfrow=c(1,1))

  plot(x=transient[1,sapply("CaD", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[1],pch=16)
  for (t in 2:9){
  points(x=transient[t,sapply("CaD", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[t],pch=16)
  }
  plot(x=transient[1,sapply("CaP", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[1],pch=16)
  for (t in 2:9){
    points(x=transient[t,sapply("CaP", grepl, colnames(transient))],
           y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[t],pch=16)
  }
  
  plot(x=transient[1,sapply("ArD", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[1],pch=16)
  for (t in 2:9){
        points(x=transient[t,sapply("ArD", grepl, colnames(transient))],
           y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[t],pch=16)
  }
  plot(x=transient[1,sapply("pH", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[1],pch=16,xlim=c(6.5,8.5))
  for (t in 2:9){
    points(x=transient[t,sapply("pH", grepl, colnames(transient))],
           y=PL$grid$x.mid,ylim=c(20,0),col=magma(9)[t],pch=16)
  }
  
  plot(x=transient[1,sapply("Rmin.f", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(2,0),col=magma(9)[1],pch=16)
  plot(x=transient[1,sapply("Rmin.s", grepl, colnames(transient))],
       y=PL$grid$x.mid,ylim=c(10,0),col=magma(9)[1],pch=16)
  
  
  print(paste("water depth=",PL$water.depth,"m"))
  print(paste("AT flux=",-(ModelResults[[i]]$output$F.TA.up+ModelResults[[i]]$output$TA.irr.int)*1e-3*1e4/365.25,"meq m-2 d-1"))
  
  #points(x=PL$water.depth,
  #       y=-(ModelResults[[i]]$output$F.TA.up+ModelResults[[i]]$output$TA.irr.int-2.*(ModelResults[[i]]$output$F.H2S.up+ModelResults[[i]]$output$H2S.irr.int))*1e-3*1e4/365.25,pch=15,col="red")
  
#}

#--------------------------------------------------
# Coastal Sand
#--------------------------------------------------

i <- 1
  
PL <- set.parameter.list(PL=NULL)
  
PL$water.depth <- water.depth[i] 
PL$irr.0       <- 1.5*365.25*(5.2 * 10^(0.76241122-0.00039724*PL$water.depth))/30.# assuming irrigation of 1 d-1 is highest possible and irrigation scales with bioturbation
PL$x.irr       <- 3.
PL$Db.0        <- 5.#5.2 * 10^(0.76241122-0.00039724*PL$water.depth) 
PL$x.L         <- 2.#1.+9.*(1.-exp(-PL$Db.0/3.))
PL$v.sed       <- 3.3 * 10^(-0.87478367-0.00043512*PL$water.depth) 
PL$Cflux       <- 5.*365.25*1e3*1e-4#4/5*40.61*PL$water.depth^(-0.571)*1e6*1e-4 # umol cm-2 yr-1
#PL$F.CaCO3     <- 2./3.*0.06*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
#PL$F.CaCO3.arg <- 1./3.*0.06*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
PL$F.CaCO3     <- 2./3.*0.5*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
PL$F.CaCO3.arg <- 1./3.*0.5*PL$Cflux # Take fixed POC:PIC ratio, based on Krumins
PL$NO3.ow      <- 0.003734 + 0.0000584*PL$water.depth
PL$CO2.ow      <- DIC.BW[i]
PL$TA.ow       <- TA.BW[i]
PL$TC          <- 15.69 - 0.054*PL$water.depth
if(PL$TC<2.) PL$TC <- 2.

# set parameters for sandy site
PL$k.f <- 10.  # yr-1 10x higher because i don't think 1 is high enough
PL$k.s <- 0.1 # yr-1
PL$p.f <- 0.8
PL$p.s <- 0.2
PL$p.v <- 0.
PL$p_fe <- 0.0

PL$irr.H2S <- 0.6
PL$irr.CO2 <- 0.2
PL$irr.NH4 <- 0.3

#PL$k.CaP <- 0.

PL$por.0   <- 0.4
PL$por.inf <- 0.4

# Update parameter list 
PL     <- set.parameter.list(PL)
SV.ini <- set.initial.state(PL)
  
PL$Run.name <- paste("Baseline_CoastalSand",PL$water.depth,sep="_")

out.temp <- steady.1D(y=SV.ini, func=DiagC.model, parms=PL, nspec=PL$N.var, pos=TRUE)
print(paste("water depth=",PL$water.depth,"m"))
print(paste("AT flux=",-((out.temp$F.TA.up+out.temp$TA.irr.int) - (out.temp$F.H2S.up+out.temp$H2S.irr.int))*1e-3*1e4/365.25,"meq m-2 d-1"))
print(paste("DIC flux=",-((out.temp$F.CO2.up+out.temp$CO2.irr.int))*1e-3*1e4/365.25,"mmol m-2 d-1"))
print(paste("O2 flux=",((out.temp$F.O2.up+out.temp$O2.irr.int) - (out.temp$F.H2S.up+out.temp$H2S.irr.int))*1e-3*1e4/365.25,"mmol m-2 d-1"))
print(Depth.Matrix[,"Coastal sand"])

x11()
plotProfile(model.natural=out.temp,
            model.chronic=NULL,
            timepoint.chronic=NULL,
            PL=PL,
            ylimz = c(20., 0),
            colors=magma(5) )
  
ModelSandResults <- list(PL=PL,output=out.temp)

save(ModelSandResults,file="Baseline_CoastalSand.Rdata")
