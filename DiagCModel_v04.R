#####################################################################################################
# Diagenetic carbon model to assess controls on sedimentary organic carbon breakdown
# Author: Sebastiaan van de Velde - Royal Belgian Institute of Natural Sciences
#                                   Université Libre de Bruxelles
#         Elwenn Eon              - Université Libre de Bruxelles
# contact: svandevelde@naturalsciences.be
#####################################################################################################

# 09/12/2022 - v01
# SvdV extended reaction network with relevant alkalinity reactions
# 12/12/2022 - v02
# SvdV added carbonate dissolution/precipation
# included pH calculation based on Munhoven, 2013 (using the functions from Meysman's unpublished 'pH' package)

#====================================================================================================
# Auxiliary functions
#====================================================================================================

library(ReacTran)
library(AquaEnv)

#----------------------------------------------------------------------------------------------------
# Numerical limit for reactions
#----------------------------------------------------------------------------------------------------

FSAT <- function(C,K,n) (C/K)^n/((C/K)^n + 1)

#====================================================================================================
# Model function
#====================================================================================================

DiagC.model <- function (t, state, PL) 
{
  with(as.list(c(PL)),{

      if (Cforcing == "constant") { # Upper boundary condition for CH2O: fixed flux
        flux.up.CH2O.f <- p.f*Cflux/106. # from C to OM flux
        flux.up.CH2O.s <- p.s*Cflux/106. # from C to OM flux
        flux.up.CH2O.v <- p.v*Cflux/106. # from C to OM flux
      } else { # Upper boundary condition for CH2O: seasonal flux
        b <- 0.0172
        h <- 49.11
        flux.up.CH2O.f <- (sin (b * (t*365.25 - h)) + 1.)*p.f*Cflux/106. # from C to OM flux
        flux.up.CH2O.s <- (sin (b * (t*365.25 - h)) + 1.)*p.s*Cflux/106. # from C to OM flux
        flux.up.CH2O.v <- (sin (b * (t*365.25 - h)) + 1.)*p.v*Cflux/106. # from C to OM flux
      }
    
    if (PL$TI.Benthos==TRUE) {
      Db.corr <- approx(x=Db.modif[,1],y=Db.modif[,2],xout=t,rule=2)$y
    } else { Db.corr <- 1.}
    
    # Initialization of state variables 
    CH2O.f <- state[1:N.gr]
    CH2O.s <- state[(N.gr+1):(2*N.gr)]
    CH2O.v <- state[(2*N.gr+1):(3*N.gr)]
    O2     <- state[(3*N.gr+1):(4*N.gr)]
    NO3    <- state[(4*N.gr+1):(5*N.gr)]
    SO4    <- state[(5*N.gr+1):(6*N.gr)]
    CO2    <- state[(6*N.gr+1):(7*N.gr)]
    NH4    <- state[(7*N.gr+1):(8*N.gr)]
    H2PO4  <- state[(8*N.gr+1):(9*N.gr)]
    H2S    <- state[(9*N.gr+1):(10*N.gr)]
    N2     <- state[(10*N.gr+1):(11*N.gr)]
    TA     <- state[(11*N.gr+1):(12*N.gr)]
    FeS2   <- state[(12*N.gr+1):(13*N.gr)]
    CH4    <- state[(13*N.gr+1):(14*N.gr)]
    CaCO3  <- state[(14*N.gr+1):(15*N.gr)]
    Ca     <- state[(15*N.gr+1):(16*N.gr)]
    CaCO3.arg <- state[(16*N.gr+1):(17*N.gr)]
    
    # Transport terms

      tran.CH2O.f <- tran.1D(C=CH2O.f, flux.up=flux.up.CH2O.f, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid,dx=grid)$dC
      tran.CH2O.s <- tran.1D(C=CH2O.s, flux.up=flux.up.CH2O.s, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid,dx=grid)$dC
      tran.CH2O.v <- tran.1D(C=CH2O.v, flux.up=flux.up.CH2O.v, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid,dx=grid)$dC
      tran.FeS2   <- tran.1D(C=FeS2,   flux.up=F.FeS2,         v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid,dx=grid)$dC
      tran.CaCO3  <- tran.1D(C=CaCO3,  flux.up=F.CaCO3,        v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid,dx=grid)$dC
      tran.CaCO3.arg  <- tran.1D(C=CaCO3.arg, flux.up=F.CaCO3.arg, v=v.grid, D=Db.grid$int*Db.corr,  VF=svf.grid,dx=grid)$dC
      
      if (is.na(O2.ds))     O2.C.down      <- O2[N.gr]     else O2.C.down     <- O2.ds 
      if (is.na(NO3.ds))    NO3.C.down     <- NO3[N.gr]    else NO3.C.down    <- NO3.ds 
      if (is.na(SO4.ds))    SO4.C.down     <- SO4[N.gr]    else SO4.C.down    <- SO4.ds  
      if (is.na(CO2.ds))    CO2.C.down     <- CO2[N.gr]    else CO2.C.down    <- CO2.ds 
      if (is.na(NH4.ds))    NH4.C.down     <- NH4[N.gr]    else NH4.C.down    <- NH4.ds 
      if (is.na(H2PO4.ds))  H2PO4.C.down   <- H2PO4[N.gr]  else H2PO4.C.down  <- H2PO4.ds 
      if (is.na(H2S.ds))    H2S.C.down     <- H2S[N.gr]    else H2S.C.down    <- H2S.ds  
      if (is.na(N2.ds))     N2.C.down      <- N2[N.gr]     else N2.C.down     <- N2.ds  
      if (is.na(TA.ds))     TA.C.down      <- TA[N.gr]     else TA.C.down     <- TA.ds  
      if (is.na(CH4.ds))    CH4.C.down     <- CH4[N.gr]    else CH4.C.down    <- CH4.ds  
      if (is.na(Ca.ds))     Ca.C.down      <- Ca[N.gr]     else Ca.C.down     <- Ca.ds  
      
      tran.O2     <- tran.1D(C=O2,   C.up=O2.ow,   C.down=O2.C.down,   D=D.O2.grid,   v=u.grid,VF=por.grid,dx=grid)$dC
      tran.NO3    <- tran.1D(C=NO3,  C.up=NO3.ow,  C.down=NO3.C.down,  D=D.NO3.grid,  v=u.grid,VF=por.grid,dx=grid)$dC
      tran.SO4    <- tran.1D(C=SO4,  C.up=SO4.ow,  C.down=SO4.C.down,  D=D.SO4.grid,  v=u.grid,VF=por.grid,dx=grid)$dC
      tran.CO2    <- tran.1D(C=CO2,  C.up=CO2.ow,  C.down=CO2.C.down,  D=D.CO2.grid,  v=u.grid,VF=por.grid,dx=grid)$dC
      tran.NH4    <- tran.1D(C=NH4,  C.up=NH4.ow,  C.down=NH4.C.down,  D=D.NH4.grid,  v=u.grid,VF=por.grid,dx=grid)$dC
      tran.H2PO4  <- tran.1D(C=H2PO4,C.up=H2PO4.ow,C.down=H2PO4.C.down,D=D.H2PO4.grid,v=u.grid,VF=por.grid,dx=grid)$dC
      tran.H2S    <- tran.1D(C=H2S,  C.up=H2S.ow,  C.down=H2S.C.down,  D=D.H2S.grid,  v=u.grid,VF=por.grid,dx=grid)$dC
      tran.N2     <- tran.1D(C=N2,   C.up=N2.ow,   C.down=N2.C.down,   D=D.N2.grid,   v=u.grid,VF=por.grid,dx=grid)$dC
      tran.TA     <- tran.1D(C=TA,   C.up=TA.ow,   C.down=TA.C.down,   D=D.TA.grid,   v=u.grid,VF=por.grid,dx=grid)$dC
      tran.CH4    <- tran.1D(C=CH4,  C.up=CH4.ow,  C.down=CH4.C.down,  D=D.CH4.grid,  v=u.grid,VF=por.grid,dx=grid)$dC
      tran.Ca     <- tran.1D(C=Ca,   C.up=Ca.ow,   C.down=Ca.C.down,   D=D.Ca.grid,   v=u.grid,VF=por.grid,dx=grid)$dC
      
    
    # Derived variables
      
      eq.obj <- pH.equilibration(S=S,T=TC,TA=TA*1e-3,SumCO2=CO2*1e-3)
      Ohm.Ca <- (Ca*1e-3*eq.obj$CO3)/Ksp.Cal
      Ohm.Ar <- (Ca*1e-3*eq.obj$CO3)/Ksp.Ara
      pH     <- eq.obj$pH
      
      # OM reaction rate

      # Constant reactivity (changes only depending on site & OC fraction)
      if (model == "m1") { 
        k.s.O2 <- k.s.base
        k.v.O2 <- k.v.base
        
        k.s <- k.s.base
        k.v <- k.v.base
      }

      # Only for the slow and very slow OC fractions: reactivity increased for aerobic respiration 
      if (model == "m2") {
        k.s.O2 <- k.O2.acc*k.s.base
        k.v.O2 <- k.O2.acc*k.v.base        
      
        k.s <- k.s.base
        k.v <- k.v.base
      }

      # Only for the slow and very slow OC fractions: reactivity dependent on the fast OC fraction concentration 
      if (model == "m3") {
        k.s <- k.s.base + (UL.SP - 1.)*k.s.base * (1. - exp(-SP.SP*CH2O.f))

        k.v <- k.v.base + (UL.SP - 1.)*k.v.base * (1. - exp(-SP.SP*CH2O.f))

        k.s.O2 <- k.s
        k.v.O2 <- k.v 
      }    

    # Reaction rates
      
      # Aerobic respiration
      AR.f <- svf.grid$mid * k.f    * CH2O.f * (O2/(Ks.O2+O2)) * (O2>0.) * FSAT(CH2O.f,C.lim,5)
      AR.s <- svf.grid$mid * k.s.O2 * CH2O.s * (O2/(Ks.O2+O2)) * (O2>0.) * FSAT(CH2O.s,C.lim,5)
      AR.v <- svf.grid$mid * k.v.O2 * CH2O.v * (O2/(Ks.O2+O2)) * (O2>0.) * FSAT(CH2O.v,C.lim,5)
      AR   <- AR.f + AR.s + AR.v

      # Nitrate reduction
      DN.f <- svf.grid$mid * k.f * CH2O.f * (Ks.O2/(Ks.O2+O2)) * (NO3/(Ks.NO3+NO3)) * (NO3>0.) * FSAT(CH2O.f,C.lim,5)
      DN.s <- svf.grid$mid * k.s * CH2O.s * (Ks.O2/(Ks.O2+O2)) * (NO3/(Ks.NO3+NO3)) * (NO3>0.) * FSAT(CH2O.s,C.lim,5)
      DN.v <- svf.grid$mid * k.v * CH2O.v * (Ks.O2/(Ks.O2+O2)) * (NO3/(Ks.NO3+NO3)) * (NO3>0.) * FSAT(CH2O.v,C.lim,5)
      DN   <- DN.f + DN.s + DN.v
      
      # Sulfate reduction
      SR.f <- svf.grid$mid * k.f * CH2O.f * (Ks.O2/(Ks.O2+O2)) * (Ks.NO3/(Ks.NO3+NO3)) * (SO4/(Ks.SO4+SO4)) * (SO4>0.) * FSAT(CH2O.f,C.lim,5)
      SR.s <- svf.grid$mid * k.s * CH2O.s * (Ks.O2/(Ks.O2+O2)) * (Ks.NO3/(Ks.NO3+NO3)) * (SO4/(Ks.SO4+SO4)) * (SO4>0.) * FSAT(CH2O.s,C.lim,5)
      SR.v <- svf.grid$mid * k.v * CH2O.v * (Ks.O2/(Ks.O2+O2)) * (Ks.NO3/(Ks.NO3+NO3)) * (SO4/(Ks.SO4+SO4)) * (SO4>0.) * FSAT(CH2O.v,C.lim,5)
      SR   <- SR.f + SR.s + SR.v

      # Methanogenesis
      MG.f <- svf.grid$mid * k.f * CH2O.f * (Ks.O2/(Ks.O2+O2)) * (Ks.NO3/(Ks.NO3+NO3)) * (Ks.SO4/(Ks.SO4+SO4)) * FSAT(CH2O.f,C.lim,5)
      MG.s <- svf.grid$mid * k.s * CH2O.s * (Ks.O2/(Ks.O2+O2)) * (Ks.NO3/(Ks.NO3+NO3)) * (Ks.SO4/(Ks.SO4+SO4)) * FSAT(CH2O.s,C.lim,5)
      MG.v <- svf.grid$mid * k.v * CH2O.v * (Ks.O2/(Ks.O2+O2)) * (Ks.NO3/(Ks.NO3+NO3)) * (Ks.SO4/(Ks.SO4+SO4)) * FSAT(CH2O.v,C.lim,5)     
      MG   <- MG.f + MG.s + MG.v

      # Mineralization rate
      Cmin.f <- AR.f + DN.f + SR.f + MG.f
      Cmin.s <- AR.s + DN.s + SR.s + MG.s  
      Cmin.v <- AR.v + DN.v + AR.v + MG.v    
      
      # Secondary redox reactions

      Nit  <- por.grid$mid * k.Nit    * O2        * NH4 * (O2>0)    * (NH4>0)  * FSAT(O2,C.lim,5)   * FSAT(NH4,C.lim,5)
      CSO  <- por.grid$mid * k.CSO    * H2S       * O2  * (H2S>0.)  * (O2>0.)  * FSAT(O2,C.lim,5)   * FSAT(H2S,C.lim,5)
      AeOM <- por.grid$mid * k.AeOM   * CH4       * O2  * (CH4>0.)  * (O2>0.)  * FSAT(CH4,C.lim,5)  * FSAT(O2,C.lim,5)     
      AOM  <- por.grid$mid * k.AOM    * CH4       * SO4 * (CH4>0.)  * (SO4>0.) * FSAT(SO4,C.lim,5)  * FSAT(CH4,C.lim,5)
      PyO  <- svf.grid$mid * k.FeS2   * FeS2      * O2  * (FeS2>0.) * (O2>0.)  * FSAT(FeS2,C.lim,5) * FSAT(O2,C.lim,5)
      
      # Calcium carbonate dissolution - precipitation
      
      CaP  <- por.grid$mid * k.CaP * abs(Ohm.Ca - 1.)^n.CaP * (Ohm.Ca > 1.0) * (Ca>0.) * (CO2>0.) * FSAT(Ca,C.lim,5) * FSAT(CO2,C.lim,5)
      
      CaD  <- svf.grid$mid * (k.CaD[2] * abs(1. - Ohm.Ca)^n.CaD[2] * (Ohm.Ca > OhmCa.crit) +
                              k.CaD[1] * abs(1. - Ohm.Ca)^n.CaD[1] * (Ohm.Ca <= OhmCa.crit)) * CaCO3 * (Ohm.Ca < 1.0) * FSAT(CaCO3,C.lim,5)       
      
      ArD  <- svf.grid$mid * (k.ArD[2] * abs(1. - Ohm.Ar)^n.ArD[2] * (Ohm.Ar > OhmAr.crit) +
                              k.ArD[1] * abs(1. - Ohm.Ar)^n.ArD[1] * (Ohm.Ar <= OhmAr.crit))* CaCO3.arg *(Ohm.Ar < 1.0) * FSAT(CaCO3.arg,C.lim,5)       
      
      # Reaction terms
      reac.CH2O.f <- (-Cmin.f) / svf.grid$mid
      reac.CH2O.s <- (-Cmin.s) / svf.grid$mid 
      reac.CH2O.v <- (-Cmin.v) / svf.grid$mid
      
      reac.O2     <- (-106*AR   - 2.*Nit              - 2.*CSO  - 2.*AeOM - 15/4*PyO)/por.grid$mid
      reac.NO3    <- (-424/5*DN + Nit                                               )/por.grid$mid
      reac.SO4    <- (-(848-742*p_fe)/(16-15*p_fe)*SR + CSO + 2.*PyO - AOM          )/por.grid$mid
      
      reac.CO2    <- (106*AR + 106*DN + 106*SR + 58*MG + AeOM + AOM - CaP + CaD + ArD)/por.grid$mid
      reac.NH4    <- (+16* (AR + DN + SR + MG) - Nit)/por.grid$mid
      reac.H2PO4  <- +1*(AR + DN + SR + MG)/por.grid$mid
      
      reac.CH4    <- (+58*MG - AeOM - AOM)/por.grid$mid
      reac.H2S    <- (+(848-1590*p_fe)/(16-15*p_fe)*SR            - CSO + AOM)/por.grid$mid
      reac.N2     <- (+212/5*DN)/por.grid$mid
      reac.TA     <- (+15*AR + 499/5*DN + (1936-1709*p_fe)/(16-15*p_fe)*SR + 15.*MG - 2.*Nit - 2*CSO - 4.*PyO + 2.*AOM - 2.*CaP + 2.*CaD + 2.*ArD)/por.grid$mid
      
      reac.FeS2   <- (+(424*p_fe)/(16-15*p_fe)*SR - PyO)/svf.grid$mid
      
      reac.CaCO3      <- ( CaP - CaD        )/svf.grid$mid
      reac.CaCO3.arg  <- (            - ArD )/svf.grid$mid
      reac.Ca         <- (-CaP + CaD  + ArD )/por.grid$mid
      
      # irrigation
      irr.O2    <- irr.grid$mid*Db.corr*irr.O2*(O2.ow - O2) #* FSAT((O2.ow - O2),C.lim,5) 
      irr.NO3   <- irr.grid$mid*Db.corr*irr.NO3*(NO3.ow - NO3)# * FSAT((NO3.ow - NO3),C.lim,5)
      irr.SO4   <- irr.grid$mid*Db.corr*irr.SO4*(SO4.ow - SO4) #* FSAT((SO4.ow - SO4),C.lim,5)
      irr.CO2   <- irr.grid$mid*Db.corr*irr.CO2*(CO2.ow - CO2) #* FSAT((CO2.ow - CO2),C.lim,5)
      irr.NH4   <- irr.grid$mid*Db.corr*irr.NH4*(NH4.ow - NH4) #* FSAT((NH4.ow - NH4),C.lim,5)
      irr.H2PO4 <- irr.grid$mid*Db.corr*irr.H2PO4*(H2PO4.ow - H2PO4) #* FSAT((H2PO4.ow - H2PO4),C.lim,5)
      irr.CH4   <- irr.grid$mid*Db.corr*irr.CH4*(CH4.ow - CH4) #* FSAT((CH4.ow - CH4),C.lim,5)
      irr.H2S   <- irr.grid$mid*Db.corr*irr.H2S*(H2S.ow - H2S) #* FSAT((H2S.ow - H2S),C.lim,5)
      irr.N2    <- irr.grid$mid*Db.corr*irr.N2*(N2.ow - N2) #* FSAT((N2.ow - N2),C.lim,5)
      irr.TA    <- irr.grid$mid*Db.corr*irr.TA*(TA.ow - TA) #* FSAT((TA.ow - TA),C.lim,5)
      irr.Ca    <- irr.grid$mid*Db.corr*irr.Ca*(Ca.ow - Ca) #* FSAT((Ca.ow - Ca),C.lim,5)
      
      # Assemble differential equations
      ddt.CH2O.f <- tran.CH2O.f + reac.CH2O.f
      ddt.CH2O.s <- tran.CH2O.s + reac.CH2O.s
      ddt.CH2O.v <- tran.CH2O.v + reac.CH2O.v
      ddt.FeS2   <- tran.FeS2   + reac.FeS2
      ddt.CaCO3  <- tran.CaCO3  + reac.CaCO3
      ddt.CaCO3.arg  <- tran.CaCO3.arg  + reac.CaCO3.arg
      
      ddt.O2     <- tran.O2     + reac.O2    + irr.O2
      ddt.NO3    <- tran.NO3    + reac.NO3   + irr.NO3
      ddt.SO4    <- tran.SO4    + reac.SO4   + irr.SO4
      
      ddt.CO2    <- tran.CO2    + reac.CO2   + irr.CO2
      ddt.NH4    <- tran.NH4    + reac.NH4   + irr.NH4
      ddt.H2PO4  <- tran.H2PO4  + reac.H2PO4 + irr.H2PO4
      
      ddt.CH4    <- tran.CH4    + reac.CH4   + irr.CH4
      ddt.H2S    <- tran.H2S    + reac.H2S   + irr.H2S
      ddt.N2     <- tran.N2     + reac.N2    + irr.N2
      ddt.TA     <- tran.TA     + reac.TA    + irr.TA
      ddt.Ca     <- tran.Ca     + reac.Ca    + irr.Ca
      
    # Fluxes in and out of the domain

      F.Cmin.f.up   <-  tran.1D(C=CH2O.f, flux.up=flux.up.CH2O.f, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.up
      F.Cmin.f.down <-  tran.1D(C=CH2O.f, flux.up=flux.up.CH2O.f, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.down
      
      F.Cmin.s.up   <-  tran.1D(C=CH2O.s, flux.up=flux.up.CH2O.s, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.up
      F.Cmin.s.down <-  tran.1D(C=CH2O.s, flux.up=flux.up.CH2O.s, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.down   

      F.Cmin.v.up   <-  tran.1D(C=CH2O.v, flux.up=flux.up.CH2O.v, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.up
      F.Cmin.v.down <-  tran.1D(C=CH2O.v, flux.up=flux.up.CH2O.v, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.down

      F.FeS2.up     <-  tran.1D(C=FeS2, flux.up=F.FeS2, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.up
      F.FeS2.down   <-  tran.1D(C=FeS2, flux.up=F.FeS2, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.down
      
      F.CaCO3.up    <-  tran.1D(C=CaCO3, flux.up=F.CaCO3, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.up
      F.CaCO3.down  <-  tran.1D(C=CaCO3, flux.up=F.CaCO3, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.down
      
      F.CaCO3.arg.up    <-  tran.1D(C=CaCO3.arg, flux.up=F.CaCO3.arg, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.up
      F.CaCO3.arg.down  <-  tran.1D(C=CaCO3.arg, flux.up=F.CaCO3.arg, v=v.grid, D=Db.grid$int*Db.corr,    VF=svf.grid, dx=grid)$flux.down
      
      F.O2.up       <-  tran.1D(C=O2,     C.up=O2.ow,             v=v.grid, D=D.O2.grid,  VF=por.grid, dx=grid)$flux.up
      F.O2.down     <-  tran.1D(C=O2,     C.up=O2.ow,             v=v.grid, D=D.O2.grid,  VF=por.grid, dx=grid)$flux.down 
      O2.irr.int    <-  sum(grid$dx*por.grid$mid*irr.O2)
      
      F.NO3.up      <-  tran.1D(C=NO3,    C.up=NO3.ow,            v=v.grid, D=D.NO3.grid,  VF=por.grid, dx=grid)$flux.up
      F.NO3.down    <-  tran.1D(C=NO3,    C.up=NO3.ow,            v=v.grid, D=D.NO3.grid,  VF=por.grid, dx=grid)$flux.down 
      NO3.irr.int   <-  sum(grid$dx*por.grid$mid*irr.NO3)
      
      F.SO4.up      <-  tran.1D(C=SO4,    C.up=SO4.ow,            v=v.grid, D=D.SO4.grid, VF=por.grid, dx=grid)$flux.up
      F.SO4.down    <-  tran.1D(C=SO4,    C.up=SO4.ow,            v=v.grid, D=D.SO4.grid, VF=por.grid, dx=grid)$flux.down 
      SO4.irr.int   <-  sum(grid$dx*por.grid$mid*irr.SO4)
      
      F.CO2.up      <-  tran.1D(C=CO2,    C.up=CO2.ow,            v=v.grid, D=D.CO2.grid, VF=por.grid, dx=grid)$flux.up
      F.CO2.down    <-  tran.1D(C=CO2,    C.up=CO2.ow,            v=v.grid, D=D.CO2.grid, VF=por.grid, dx=grid)$flux.down 
      CO2.irr.int   <-  sum(grid$dx*por.grid$mid*irr.CO2)
      
      F.NH4.up      <-  tran.1D(C=NH4,    C.up=NH4.ow,            v=v.grid, D=D.NH4.grid, VF=por.grid, dx=grid)$flux.up
      F.NH4.down    <-  tran.1D(C=NH4,    C.up=NH4.ow,            v=v.grid, D=D.NH4.grid, VF=por.grid, dx=grid)$flux.down 
      NH4.irr.int   <-  sum(grid$dx*por.grid$mid*irr.NH4)
      
      F.H2PO4.up    <-  tran.1D(C=H2PO4,  C.up=H2PO4.ow,          v=v.grid, D=D.H2PO4.grid, VF=por.grid, dx=grid)$flux.up
      F.H2PO4.down  <-  tran.1D(C=H2PO4,  C.up=H2PO4.ow,          v=v.grid, D=D.H2PO4.grid, VF=por.grid, dx=grid)$flux.down 
      H2PO4.irr.int <-  sum(grid$dx*por.grid$mid*irr.H2PO4)
      
      F.CH4.up      <-  tran.1D(C=CH4,    C.up=CH4.ow,            v=v.grid, D=D.CH4.grid, VF=por.grid, dx=grid)$flux.up
      F.CH4.down    <-  tran.1D(C=CH4,    C.up=CH4.ow,            v=v.grid, D=D.CH4.grid, VF=por.grid, dx=grid)$flux.down 
      CH4.irr.int   <-  sum(grid$dx*por.grid$mid*irr.CH4)
      
      F.H2S.up      <-  tran.1D(C=H2S,    C.up=H2S.ow,            v=v.grid, D=D.H2S.grid, VF=por.grid, dx=grid)$flux.up
      F.H2S.down    <-  tran.1D(C=H2S,    C.up=H2S.ow,            v=v.grid, D=D.H2S.grid, VF=por.grid, dx=grid)$flux.down 
      H2S.irr.int   <-  sum(grid$dx*por.grid$mid*irr.H2S)
      
      F.N2.up       <-  tran.1D(C=N2,     C.up=N2.ow,             v=v.grid, D=D.N2.grid, VF=por.grid, dx=grid)$flux.up
      F.N2.down     <-  tran.1D(C=N2,     C.up=N2.ow,             v=v.grid, D=D.N2.grid, VF=por.grid, dx=grid)$flux.down 
      N2.irr.int    <-  sum(grid$dx*por.grid$mid*irr.N2)
      
      F.TA.up       <-  tran.1D(C=TA,     C.up=TA.ow,             v=v.grid, D=D.TA.grid, VF=por.grid, dx=grid)$flux.up
      F.TA.down     <-  tran.1D(C=TA,     C.up=TA.ow,             v=v.grid, D=D.TA.grid, VF=por.grid, dx=grid)$flux.down 
      TA.irr.int    <-  sum(grid$dx*por.grid$mid*irr.TA)
      
      F.Ca.up       <-  tran.1D(C=Ca,     C.up=Ca.ow,             v=v.grid, D=D.Ca.grid, VF=por.grid, dx=grid)$flux.up
      F.Ca.down     <-  tran.1D(C=Ca,     C.up=Ca.ow,             v=v.grid, D=D.Ca.grid, VF=por.grid, dx=grid)$flux.down 
      Ca.irr.int   <-  sum(grid$dx*por.grid$mid*irr.Ca)
      
    # Return the total rates of change, fluxes and reaction rates

      return(list(
        # Total rates of change
        c(ddt.CH2O.f, ddt.CH2O.s, ddt.CH2O.v, ddt.O2, ddt.NO3, ddt.SO4, ddt.CO2, ddt.NH4, ddt.H2PO4, ddt.H2S, 
          ddt.N2, ddt.TA, ddt.FeS2, ddt.CH4, ddt.CaCO3, ddt.Ca, ddt.CaCO3.arg),
       
        "pH"=pH,

        # Mineralization
        "Rmin.f" = Cmin.f, "Rmin.s" = Cmin.s, "Rmin.v" = Cmin.v,

        # Reaction rates
        "AR"  = AR,  "DN" = DN,   "SR"  = SR,    "MG"  = MG, 
        "CSO" = CSO, "AOM" = AOM, "AeOM" = AeOM, "Nit" = Nit,
        "PyO" = PyO, 'CaP' = CaP, "CaD" = CaD,   "ArD" = ArD,

        # Fluxes
        F.Cmin.f.up   = F.Cmin.f.up,   F.Cmin.s.up   = F.Cmin.s.up,   F.Cmin.v.up   = F.Cmin.v.up,   F.FeS2.up   = F.FeS2.up,   F.CaCO3.up   = F.CaCO3.up,   F.CaCO3.arg.up   = F.CaCO3.arg.up,
        F.Cmin.f.down = F.Cmin.f.down, F.Cmin.s.down = F.Cmin.s.down, F.Cmin.v.down = F.Cmin.v.down, F.FeS2.down = F.FeS2.down, F.CaCO3.down = F.CaCO3.down, F.CaCO3.arg.down = F.CaCO3.arg.down,
        
        F.O2.up   = F.O2.up,   F.NO3.up   = F.NO3.up,   F.SO4.up   = F.SO4.up,   F.CO2.up   = F.CO2.up,   F.NH4.up   = F.NH4.up,   F.H2PO4.up   = F.H2PO4.up, 
        F.O2.down = F.O2.down, F.NO3.down = F.NO3.down, F.SO4.down = F.SO4.down, F.CO2.down = F.CO2.down, F.NH4.down = F.NH4.down, F.H2PO4.down = F.H2PO4.down,
        
        F.CH4.up   = F.CH4.up,   F.H2S.up   = F.H2S.up,   F.N2.up   = F.N2.up,   F.TA.up   = F.TA.up,   F.Ca.up   = F.Ca.up, 
        F.CH4.down = F.CH4.down, F.H2S.down = F.H2S.down, F.N2.down = F.N2.down, F.TA.down = F.TA.down, F.Ca.down = F.Ca.down,
        
        O2.irr.int  = O2.irr.int,  NO3.irr.int = NO3.irr.int, SO4.irr.int = SO4.irr.int, CO2.irr.int = CO2.irr.int, NH4.irr.int = NH4.irr.int, H2PO4.irr.int = H2PO4.irr.int,
        CH4.irr.int = CH4.irr.int, H2S.irr.int = H2S.irr.int, N2.irr.int  = N2.irr.int,  TA.irr.int  = TA.irr.int,  Ca.irr.int  = Ca.irr.int
        
        )
      )
    
  })}
