###############################################################################
# Auxilliary functions for pH calculations
###############################################################################

#require(AquaEnv)

# The file "pH_estimation.f90" must be compiled to produce the file
# "pH_estimation.dll" 
#system("RCMD SHLIB pH_estimation.f90")

# The file "pH_estimation.dll" must be dynamically loaded before running this
# script
#dyn.load("pH_estimation.dll")
#dyn.load("../package/pH_estimation.dll")
dyn.load("../Functions/pH.dll")

#=============================================================================
# Full equilibration (pH calculation + acid-base speciation) 
#=============================================================================

pH.equilibration <- function(S=35, T=25,
K_W=NULL, K_CO2=NULL, K_HCO3=NULL, K_BOH3=NULL, K_NH4=NULL, K_H3PO4=NULL, K_H2PO4=NULL, K_HPO4=NULL, K_SiOH4=NULL, K_H2S=NULL, K_HSO4=NULL, K_HF=NULL,
TA, pH=rep(7,length(SumCO2)), SumCO2, SumBOH3=NULL, SumNH4=rep(0,length(SumCO2)), SumH3PO4=rep(0,length(SumCO2)), SumSiOH4=rep(0,length(SumCO2)), SumH2S=rep(0,length(SumCO2)), SumH2SO4=NULL, SumHF=NULL,
Precision = 1E-8, Href = 1)
{

  # Make sure the function arguments have right type
  
  N <- as.integer(length(SumCO2))
  Iterations <- as.integer(rep(0,N))
  Precision <- as.double(Precision)
  Href <- as.double(Href)
  
  # Initialisation of variables 
  
  H <- Href*10^(-pH)
  
  OH <- CO2 <- HCO3 <- CO3 <- BOH3 <- BOH4 <- NH3 <- NH4 <- rep(0,N)
  H3PO4 <- H2PO4 <- HPO4 <- PO4 <- SiOH4 <- SiOOH3 <- rep(0,N)
  H2S <- HS <- HSO4 <- SO4 <- HF <- F <- rep(0,N)
  
  # Concentrations of total species 
  # Dickson et al. [2007, chap. 5, p. 10] 
  
  if (is.null(SumH2SO4)) SumH2SO4 <- rep(ConcRelCl$SO4/MeanMolecularMass$SO4*(S/1.80655),N)
  if (is.null(SumBOH3)) SumBOH3 <- rep(ConcRelCl$B/MeanMolecularMass$B*(S/1.80655),N)
  if (is.null(SumHF)) SumHF <- rep(ConcRelCl$F/MeanMolecularMass$F*(S/1.80655),N)
  
  # Calculate equilibrium constants via AquaEnv in reference situation
  if (is.null(K_W))     K_W      <- K_W(S,T,SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href^2
  if (is.null(K_BOH3))  K_BOH3   <- K_BOH3(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_HPO4))  K_HPO4   <- K_HPO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H2PO4)) K_H2PO4  <- K_H2PO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H3PO4)) K_H3PO4  <- K_H3PO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_SiOH4)) K_SiOH4  <- K_SiOH4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H2S))   K_H2S    <- K_H2S(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_HSO4))  K_HSO4   <- K_HSO4(S, T)[[1]]*Href
  if (is.null(K_HF))    K_HF     <- K_HF(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_NH4))   K_NH4    <- K_NH4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_CO2))   K_CO2    <- K_CO2(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href, k1k2="millero")[[1]]*Href
  if (is.null(K_HCO3))  K_HCO3   <- K_HCO3(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href, k1k2="millero")[[1]]*Href
  
# Call Fortan routine

  res <- .Fortran("pH_equilibration_vect", N, 
    K_W, K_CO2, K_HCO3, K_BOH3, K_NH4, K_H3PO4, K_H2PO4, K_HPO4,
    K_SiOH4, K_H2S, K_HSO4, K_HF,
    TA, SumCO2, SumBOH3,
    SumNH4, SumH3PO4, SumSiOH4, SumH2S, SumH2SO4, SumHF,
    Href, Iterations, Precision,
    pH, H, OH, CO2, HCO3, CO3, BOH3, BOH4, NH4, NH3,                       
    H3PO4, H2PO4, HPO4, PO4, SiOH4, SiOOH3, H2S, HS, HSO4, SO4, HF, F)
  
names(res) <- c("N","K_W","K_CO2","K_HCO3","K_BOH3","K_NH4",
"K_H3PO4","K_H2PO4","K_HPO4","K_SiOH4",
"K_H2S","K_HSO4","K_HF",
"TA","SumCO2","SumBOH3",
"SumNH4","SumH3PO4","SumSiOH4","SumH2S","SumH2SO4","SumHF",
"Href","Iterations","Precision",
"pH", "H", "OH", "CO2", "HCO3", "CO3", "BOH3", "BOH4", "NH4", "NH3",                       
"H3PO4", "H2PO4", "HPO4", "PO4", "SiOH4", "SiOOH3", "H2S", "HS", "HSO4", "SO4", "HF", "F")

return(res)
}

#=============================================================================
# pH calculation 
#=============================================================================

pH.calculation <- function(S=35, T=25,
K_W=NULL, K_CO2=NULL, K_HCO3=NULL, K_BOH3=NULL, K_NH4=NULL, K_H3PO4=NULL, K_H2PO4=NULL, K_HPO4=NULL, K_SiOH4=NULL, K_H2S=NULL, K_HSO4=NULL, K_HF=NULL,
TA, pH=7, SumCO2, SumBOH3=NULL, SumNH4=0, SumH3PO4=0, SumSiOH4=0, SumH2S=0, SumH2SO4=NULL, SumHF=NULL,
Precision = 1E-8, Href = 1)
{

# Make sure the function arguments have right type

  N <- as.integer(length(SumCO2))
  Iterations <- as.integer(rep(0,N))
  Precision <- as.double(Precision)
  Href <- as.double(Href)
  

# Initialisation of variables 

H <- Href*10^(-pH)

# Concentrations of total species 
# Dickson et al. [2007, chap. 5, p. 10] 

if (is.null(SumH2SO4)) SumH2SO4 <- rep(ConcRelCl$SO4/MeanMolecularMass$SO4*(S/1.80655),N)*Href
if (is.null(SumBOH3)) SumBOH3 <- rep(ConcRelCl$B/MeanMolecularMass$B*(S/1.80655),N)*Href
if (is.null(SumHF)) SumHF <- rep(ConcRelCl$F/MeanMolecularMass$F*(S/1.80655),N)*Href

  # Calculate equilibrium constants via AquaEnv in reference situation
  if (is.null(K_W))     K_W      <- K_W(S,T,SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href^2
  if (is.null(K_BOH3))  K_BOH3   <- K_BOH3(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_HPO4))  K_HPO4   <- K_HPO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H2PO4)) K_H2PO4  <- K_H2PO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H3PO4)) K_H3PO4  <- K_H3PO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_SiOH4)) K_SiOH4  <- K_SiOH4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H2S))   K_H2S    <- K_H2S(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_HSO4))  K_HSO4   <- K_HSO4(S, T)[[1]]*Href
  if (is.null(K_HF))    K_HF     <- K_HF(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_NH4))   K_NH4    <- K_NH4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_CO2))   K_CO2    <- K_CO2(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href, k1k2="millero")[[1]]*Href
  if (is.null(K_HCO3))  K_HCO3   <- K_HCO3(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href, k1k2="millero")[[1]]*Href
  
  
# Call Fortran routine

res <- .Fortran("pH_calculation_vect", N, 
K_W, K_CO2, K_HCO3, K_BOH3, K_NH4, K_H3PO4, K_H2PO4, K_HPO4,
K_SiOH4, K_H2S, K_HSO4, K_HF,
TA, H, pH, SumCO2, SumBOH3,
SumNH4, SumH3PO4, SumSiOH4, SumH2S, SumH2SO4, SumHF,
Href, Iterations, Precision)

names(res) <- c("N","K_W","K_CO2","K_HCO3","K_BOH3","K_NH4",
"K_H3PO4","K_H2PO4","K_HPO4","K_SiOH4",
"K_H2S","K_HSO4","K_HF",
"TA", "H", "pH", "SumCO2","SumBOH3",
"SumNH4","SumH3PO4","SumSiOH4","SumH2S","SumH2SO4","SumHF",
"Href","Iterations","Precision")

return(res)
}

#=============================================================================
# acid-base speciation 
#=============================================================================

pH.speciation <- function(S=35, T=25,
K_W=NULL, K_CO2=NULL, K_HCO3=NULL, K_BOH3=NULL, K_NH4=NULL, K_H3PO4=NULL, K_H2PO4=NULL, K_HPO4=NULL, K_SiOH4=NULL, K_H2S=NULL, K_HSO4=NULL, K_HF=NULL,
H, SumCO2, SumBOH3=NULL, SumNH4=NULL, SumH3PO4=NULL, SumSiOH4=NULL, SumH2S=NULL, SumH2SO4=NULL, SumHF=NULL, Href=1)
{

# Make sure the function arguments have right type

N <- as.integer(length(SumCO2))

OH <- CO2 <- HCO3 <- CO3 <- BOH3 <- BOH4 <- NH3 <- NH4 <- rep(0,N)
H3PO4 <- H2PO4 <- HPO4 <- PO4 <- SiOH4 <- SiOOH3 <- rep(0,N)
H2S <- HS <- HSO4 <- SO4 <- HF <- F <- rep(0,N)

# Concentrations of total species 
# Dickson et al. [2007, chap. 5, p. 10] 

if (is.null(SumNH4)) SumNH4 <- rep(0,N)
if (is.null(SumH3PO4)) SumH3PO4 <- rep(0,N)
if (is.null(SumSiOH4)) SumSiOH4 <- rep(0,N)
if (is.null(SumH2S)) SumH2S <- rep(0,N)

if (is.null(SumH2SO4)) SumH2SO4 <- Href*rep(ConcRelCl$SO4/MeanMolecularMass$SO4*(S/1.80655),N)*Href
if (is.null(SumBOH3)) SumBOH3 <- Href*rep(ConcRelCl$B/MeanMolecularMass$B*(S/1.80655),N)*Href
if (is.null(SumHF)) SumHF <- Href*rep(ConcRelCl$F/MeanMolecularMass$F*(S/1.80655),N)*Href

# Calculate equilibrium constants via AquaEnv in reference situation

  # Calculate equilibrium constants via AquaEnv in reference situation
  if (is.null(K_W))     K_W      <- K_W(S,T,SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href^2
  if (is.null(K_BOH3))  K_BOH3   <- K_BOH3(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_HPO4))  K_HPO4   <- K_HPO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H2PO4)) K_H2PO4  <- K_H2PO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H3PO4)) K_H3PO4  <- K_H3PO4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_SiOH4)) K_SiOH4  <- K_SiOH4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_H2S))   K_H2S    <- K_H2S(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_HSO4))  K_HSO4   <- K_HSO4(S, T)[[1]]*Href
  if (is.null(K_HF))    K_HF     <- K_HF(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_NH4))   K_NH4    <- K_NH4(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href)[[1]]*Href
  if (is.null(K_CO2))   K_CO2    <- K_CO2(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href, k1k2="millero")[[1]]*Href
  if (is.null(K_HCO3))  K_HCO3   <- K_HCO3(S, T, SumH2SO4=SumH2SO4/Href, SumHF=SumHF/Href, k1k2="millero")[[1]]*Href
  
  
  
# Call to Fortran routine 

res <- .Fortran("pH_speciation_vect", N, 
K_W, K_CO2, K_HCO3, K_BOH3, K_NH4, K_H3PO4, K_H2PO4, K_HPO4,
K_SiOH4, K_H2S, K_HSO4, K_HF,
H, SumCO2, SumBOH3,
SumNH4, SumH3PO4, SumSiOH4, SumH2S, SumH2SO4, SumHF,
OH, CO2, HCO3, CO3, BOH3, BOH4, NH4, NH3,                       
H3PO4, H2PO4, HPO4, PO4, SiOH4, SiOOH3, H2S, HS, HSO4, SO4, HF, F)

names(res) <- c("N","K_W","K_CO2","K_HCO3","K_BOH3","K_NH4",
"K_H3PO4","K_H2PO4","K_HPO4","K_SiOH4",
"K_H2S","K_HSO4","K_HF",
"H", "SumCO2","SumBOH3",
"SumNH4","SumH3PO4","SumSiOH4","SumH2S","SumH2SO4","SumHF",
"OH", "CO2", "HCO3", "CO3", "BOH3", "BOH4", "NH4", "NH3",                       
"H3PO4", "H2PO4", "HPO4", "PO4", "SiOH4", "SiOOH3", "H2S", "HS", "HSO4", "SO4", "HF", "F")

return(res)
}

