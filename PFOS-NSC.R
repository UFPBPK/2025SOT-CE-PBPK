#------------------------------------------------------------------------------- 
# Build a PBPK model for PFOS model;            
# Model code written by Wei-Chun Chou on Dec. 16, 2023; reviewed and updated by Zhoumeng Lin on Jan. 10, 2024               
# Winter school hands-on exercise
# Note: Physiological and chemical-specific parameters were collected from ---
# Loccisano et al. (2012); table 1; doi: 10.1016/j.reprotox.2011.04.006
#-------------------------------------------------------------------------------

## install R packages 
library(ggplot2)      # library for plotting
library(reshape2)     # library for reshaping data (tall-narrow <-> short-wide)
library(deSolve)      # library for solving differential equations
library(minpack.lm)   # library for least squares fit using levenberg-marquart algorithm
library(FME)          # library for MCMC simulation and model fitting
library(dplyr)        # library for Grammar of Data Manipulation 
library(truncnorm)    # library for the truncated normal distribution function
library(EnvStats)     # library for Environmental Statistics, Including US EPA Guidance
library(invgamma)     # library for inverse gamma distribution function
library(foreach)      # library for parallel computing
library(doParallel)   # library for parallel computing
library(bayesplot)    # library for Plotting Bayesian models
library(pheatmap)     # library for heatmap

################################################## Define fixed parameters ##########################################################
#! Simulation (dosing) parameters
#! Oral Gavage Dosing
PDOSEoral = 15        # (mg/kg-BW), Oral dose concentration 
OGtime = 0.01         # (hrs), Length of oral gavage 
tinterval = 24        # (hrs), Varied dependent on the exposure paradigm
tdose = 1             # (no units), Dose times

#! IV Dosing
IVdoseC= 0            # (mg/ kg), IV dose concentration
IVTime = 0.01         # (hrs), length of IV infusion

#! Physiological parameters
BW  = 0.4             # (kg), Body weight; EPA Factors Handbook, 2011
MKC = 0.0084          # (% BW), Fraction mass of kidney ; Brown 1997                                           
MLC = 0.034           # (% BW), Fraction mass of Liver ; Brown 1997 

#! Cardiac Output and Blood flow
QCC = 14              # (L/h/kg^0.75), Cardiac output; Brown 1997, Forsyth 1968
QLC = 0.183           # (% QC), Fraction blood flow to liver; Brown 1997, Fisher 2000 
QKC = 0.141           # (% QC), Fraction blood flow to kidney; Brown 1997, Forsyth 1968
Htc = 0.46            # (no units), Hematocrit for the human; Davies and Morris 1993, Brown 1997

#! Tissue Volume
VPlasC =0.0312        # (L/kg BW), Fraction vol. of plasma  ; Davies 1993
VLC = 0.035           # (L/kg BW), Fraction vol. of liver   ; Brown 1997
VKC = 0.0084          # (L/kg BW), Fraction vol. of kidney  ; Brown 1997
VfilC = 8.4e-4        # (L/kg BW), Fraction vol. of filtrate
VPTCC = 1.35e-4       # (L/kg kidney), Vol. of proximal tubule cells (L/g kidney) (60 million PTC cells/gram kidney, 1 PTC = 2250 um3) CHECK

#! Chemical specific parameters
MW = 500.13                       # (g/ mol), PFOS molecular mass 
Free = 0.02                       # (no units),Free fraction in plasma (male); Free = 0.09 (female)
Vmax_baso_invitro = 393.45        # (pmol/mg protein/ min), Vmax of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 from Nakagawa, 2007
Km_baso = 27.2                    # (mg/ L), Km of basolateral transporter; averaged in vitro value of rOAT1 and rOAT3 from Nakagawa, 2007
Vmax_apical_invitro = 9300        # (pmol/mg protein/ min), Vmax of apical transporter; averaged in vitro value of Oatp1a1 from Weaver, 2010
Km_apical= 52.3                   # (mg/ L), Km of apical transporter, in vitro value for Oatp1a1 from Weaver, 2010 
RAFbaso = 4.07                    # (no units), Relative activity factor, basolateral transporter (male) (fit to data); 0.01356 (female)
RAFapi = 35                       # (no units), Relative activity factor, apical transporter (fit to data); 0.001356 (female)
protein = 2.0e-6                  # (mg protein/proximal tubuel cell)      | Amount of protein in proximal tubule cells
GFRC = 62.1                       # (L/hr/kg kidney), Glomerular filtration rate (male); 41.04 (female); Corley, 2005

#! Partition coefficients
PL = 3.72                         # (no unit), liver-to-blood
PK = 0.8                          # (no unit), Kidney-to-blood
PR = 0.2                          # (no unit), liver-to-blood

#! rate constants
Kdif = 0.001                      # (L/h), Diffusion rate from proximal tubule cells
Kabsc = 2.12                      # (1/(h*BW^0.25)), Rate of absorption of chemical from small intestine to liver (fit to data)                        
KunabsC = 7.06e-5                 # (1/(h*BW^0.25)), Rate of un-absorbed dose to appear in feces (fit to data)
GEC = 0.54                        # (1/(h*BW^0.25)), Gastric emptying time; from yang, 2013
K0C = 1                           # (1/(h*BW^0.25)), Rate of uptake from the stomach into the liver (fit to data)

KeffluxC = 2.49                   # (1/(h*BW^0.25)), Rate of clearance of PFOA from proximal tubule cells into blood
KbileC = 0.25                     # (1/(h*BW^0.25)), Biliary elimination rate (male); liver to feces storage (fit to data)  
KurineC = 0.32                    # (1/(h*BW^0.25)), Rate of urine elimination from urine storage (male) (fit to data)


## Create a function to run PBPK model
pbpkmodel <- function(Time, State, Parameters) {
  with (as.list(State), {
    with (as.list(Parameters), {
      
      ############################ Dosing equation ###############################
      #! IV Dosing
      IVon = ifelse (T >= IVTime, 1 ,1)                   # Switch to turn IV dose on. work for single dose
      IVdose = IVdoseC*BW                                 # (mg), Iv Dose
      IVR = IVon*IVdose/IVTime                            # (mg/hr), IV Dose rate
      dAIV = IVR                            
      
      #! Oral Dosing
      Odose=PDOSEoral*BW
      oralR = Odose/OGtime
      RDOSEoral <- oralR * (Time <= tdose * tinterval) * (Time %% tinterval < OGtime)
      dAOG = RDOSEoral  
      
      ########################################################## Scaled Parameters ##########################################################
      #! Cardiac output and blood flows
      QC = QCC*(BW**0.75)*(1-Htc)                            # (L/h), Cardiac output; adjusted for plasma
      QK = (QKC*QC)                                          # (L/h), Plasma flow to kidney
      QL = (QLC*QC)                                          # (L/h), Plasma flow to liver
      QR = QC-QK-QL                                          # (L/h), Balance check of blood flows; should equal zero
      QBal = QC - (QK+QL+QR)
      
      #! Tissue Volumes
      MK = MKC*BW*1000                                       # (g), Mass of the Kidney                                       
      VPTC = MK*VPTCC                                        # (L), Volume of proximal tubule cells
      VPlas = VPlasC*BW                                      # (L), Volume of plasma
      VK=VKC*BW                                              # (L), Volume of Kidney
      VKb = VK*0.16                                          # (L), Volume of blood in the Kidney; fraction blood volume of kidney (0.16) from Brown, 1997
      Vfil = VfilC * BW                                      # (L), Volume of filtrate
      VL = VLC*BW                                            # (L), Volume of liver
      VR=(0.93*BW)- VPlas - VK - VL-Vfil-VPTC                # (L), Volume of remaining tissue
      VBal = (0.93*BW)- (VR+VL+VK+VPlas+Vfil+VPTC)           # Balance check of tissue volumes; should equal zero
      ML = MLC*BW*1000                                       # (g), Mass of the liver
      
      #! Kidney parameters
      PTC = MKC*6e7*1000                                     # (cells/kg BW), Number of PTC (cells/kg BW) (based on 60 million PTC/gram kidney)
      MPTC = VPTC*1000                                       # (g), mass of the proximal tubule cells (assuming density 1 kg/L)
      Vmax_basoC=(Vmax_baso_invitro*RAFbaso*PTC*protein*60*(MW/1e12)*1000)       # (mg/h/kg BW^0.75)      | Vmax of basolateral transporters (average Oat1 and Oat3)
      Vmax_apicalC = (Vmax_apical_invitro*RAFapi*PTC*protein*60*(MW/1e12)*1000)  # (mg/h/kg BW^0.75)      | Vmax of apical transporters in in vitro studies (Oatp1a1)
      Vmax_baso=Vmax_basoC*BW**0.75                          # (mg/h)
      Vmax_apical=Vmax_apicalC*BW**0.75                      # (mg/h)
      Kbile = KbileC*BW**(-0.25)                             # (1/h), Biliary elimination; liver to feces storage
      Kurine= KurineC*BW**(-0.25)                            # (1/h), Urinary elimination; from filtrate
      Kefflux = KeffluxC*BW**(-0.25)                         # (1/h), Efflux clearance rate from PTC to blood 
      GFR = GFRC*(MK/1000)                                   # (L/h), Glomerular filtration rate, scaled to mass of kidney 
      
      #GI tract parameters
      Kabs = Kabsc*BW**(-0.25)                               # (1/h), rate of absorption of chemical from small intestine to liver
      Kunabs = KunabsC*BW**(-0.25)                           # (1/h), rate of un-absorbed dose to appear in feces
      GE = GEC*BW**(-0.25)                                   # (1/h), Gastric emptying time 
      K0 = K0C*BW**(-0.25)                                   # (1/h), Rate of uptake from the stomach into the liver
      
      
      ########################### Compartment Equations ############################
      #! Concentration of the chemical in vein compartment
      CVR = AR/(VR*PR)                                      # (mg/L), Concentration in the venous blood leaving the Rest of body
      CVK = AKb/ (VKb*PK)                                   # (mg/L), Concentration in the venous blood leaving the kidney
      CVL = AL/(VL*PL)                                      # (mg/L), Concentration in the venous blood leaving the liver 
      
      #! Concentration of the PFOA in the plasma
      CA_free = APlas_free/VPlas                            # (mg/L), Free PFOAC concentration in the plasma
      CA = CA_free/Free                                     # (mg/L), Concentration of total PFOA in plasma
      
      #! Rest of Body
      RR = QR*(CA-CVR)*Free                                 # (mg/h), Rate of change in rest of body
      dAR = RR                                              # (mg), Amount in rest o body 
      CR = AR/VR                                            # (mg/L),  Concentration of PFOA in the compartment of rest of body
      dAUCCR = CR                                           # (mg/L*hr), Area under curve of PFOA in the compartment of rest of body
      
      #! Kidney compartment plus 2 sub-compartment (Proximal Tubule cells: PTC, Filtrate: Fil)
      #! Concentration in kidney, PTC and fil
      CKb = AKb/VKb                                         # (mg/L), Concentration of PFOA in Kidney blood
      CPTC = APTC/VPTC                                      # (mg/L), Concentration of PFOA in PTC blood
      Cfil = AFil/Vfil                                      # (mg/L), Concentration of PFOA in FIL blood
      
      #! Virtual compartment: 
      #! Basolateral (baso) 
      #! transport, Diffusion (dif), 
      #! Apical (apical) transport, and 
      #! efflux clearance (efflux)
      #! clearance (CL) via glomerular filtration 
      RA_baso = (Vmax_baso*CKb)/(Km_baso+CKb)               # (mg/h), Rate of basolateral transporter 
      dA_baso = RA_baso                                     # (mg), Amount of basolateral transporter
      RA_apical = (Vmax_apical*Cfil)/(Km_apical + Cfil)     # (mg/h), Rate of apical transporter 
      dA_apical = RA_apical                                 # (mg), Amount of apical transporter
      Rdif = Kdif*(CKb - CPTC)                              # (mg/h), Rate of diffusion from into the PTC
      dAdif = Rdif                                          # (mg), Amount moved via glomarular filtration
      RAefflux = Kefflux*APTC                               # (mg/h), Rate of efflux clearance rate from PTC to blood
      dAefflux = RAefflux                                   # (mg), Amount of efflux clearance rate from PTC to blood
      RCI = CA*GFR*Free                                     # (mg/h), Rate of clearance (CL) to via glormerular filtration (GFR) 
      dACI = RCI                                            # (mg), Amount of clearance via GFR
      
      #! Proximal Tubule Cells (PTC)
      RPTC = Rdif + RA_apical + RA_baso - RAefflux          # (mg/h), Rate of change in PTC 
      dAPTC = RPTC                                          # (mg), Amount moved in PTC
      CPTC = APTC/ VPTC                                     # (mg/L), Concentration in PTC
      dAUCCPTC = CPTC                                       # (mg/L*h), Area under curve of PFOA in the compartment of PTC
      
      #! Proximal Tubule Lumen/ Filtrate (Fil)
      Rfil = CA*GFR*Free - RA_apical - AFil*Kurine          # (mg/h), Rate of change in Fil
      dAFil = Rfil                                          # (mg), Amount moved in Fil
      Cfil = AFil/Vfil                                      # (mg/L), Concentration in Fil
      dAUCfil = Cfil                                        # (mg/L*h), Area under curve of PFOA in the compartment of Fil                                    
      
      #! Kidney compartment 
      RKb = QK*(CA-CVK)*Free - CA*GFR*Free - Rdif - RA_baso # (mg/h), Rate of change in Kidney compartment
      dAKb = RKb                                            # (mg), Amount in kidney compartment
      CKb = AKb/VKb                                         # (mg/l), Concentration in Kidney compartment
      dAUCKb = CKb                                          # (mg/L*h), Area under curve of PFOA in the Kidney compartment
      
      #! Gastrointestinal (GI) tract
      #! Stomach compartment
      RST = RDOSEoral- K0*AST - GE*AST                      # (mg/h), Rate of change in Stomach compartment
      dAST = RST                                            # (mg), Amount in Stomach
      RabsST = K0*AST                                       # (mg/h), Rate of absorption in the stomach
      dAabsST = RabsST                                      # (mg),  Amount absorbed in the stomach
      
      #! Small intestine compartment
      RSI = GE*AST - Kabs*ASI - Kunabs*ASI                  # (mg/h), Rate of change in Small intestine compartment
      dASI = RSI                                            # (mg), Amount in Small intestine
      RabsSI = Kabs*ASI                                     # (mg/h), Rate of absorption in the Small intestine
      dAabsSI = RabsSI                                      # (mg), Amount absorbed in the Small intestine
      Total_oral_uptake = AabsSI + AabsST                   # (mg), Total oral uptake in the GI
      
      #! Liver compartment
      RL = QL*(CA-CVL)*Free - Kbile*AL + Kabs*ASI + K0*AST  # (mg/h), Rate of change in liver compartment
      dAL = RL                                              # (mg), Amount in liver compartment
      CL = AL/VL                                            # (mg/L), Concentration in liver compartment
      dAUCCL = CL                                           # (mg/L*h), Area under curve of PFOA in liver compartment
      
      #! Plasma
      RPlas_free = (QR*CVR*Free) + (QK*CVK*Free) +          # (mg/h), Rate of change in the plasma
        (QL*CVL*Free) - (QC*CA*Free) + IVR + RAefflux
      dAPlas_free = RPlas_free                              # (mg), Amount in the plasma
      dAUCCA_free = CA_free                                 # (mg/L*h), Area under curve of PFOA in liver compartment
      dAUCCA = CA
      
      #! Urine elimination
      Rurine = Kurine*AFil                                  # (mg/h), Rate of change in urine
      dAurine = Rurine                                      # (mg), Amount in urine
      percentOD_in_urine = (Aurine/Odose)*100               # (%), Percent of oral dose in the urine
      
      #! Biliary excretion
      Abile = Kbile*AL                                      # (mg), Amount of PFOA in bile excretion
      amount_per_gram_liver = (AL/ML)*1000                  # (ug/g), Amount of PFOA in liver per gram liver
      
      #! Feces compartment
      Rfeces = Kbile*AL + Kunabs*ASI                        # (mg/h), Rate of change in feces compartment
      dAfeces = Rfeces                                      # (mg), Amount of the feces compartment
      percentOD_in_feces = (Afeces/Odose)*100               # (%), Percent of the oral dose in the feces
      
      
      
      #! Mass Balance Check
      Atissue = APlas_free + AR + AKb + AFil + APTC + AL + AST + ASI
      Aloss = Aurine + Afeces
      Atotal = Atissue + Aloss
      BAL = AOG - Atotal
      
      
      list(c(dAIV,dAOG, dAR, dAKb, dACI, dAdif, dA_baso, dA_apical, 
             dAefflux, dAPTC, dAFil, dAurine, dAST, dAabsST, 
             dASI, dAabsSI, dAfeces, dAL, dAPlas_free,
             dAUCCA_free, dAUCCA,dAUCCR, dAUCCPTC, dAUCfil, dAUCKb,
             dAUCCL), 
           c("CA_free"=CA_free,"CA"=CA))
      
    }) # end with parameters
  }) # end with y
} # end bd.model


State <- c(AIV=0,AOG=0, AR=0, AKb=0, ACI=0, Adif=0, A_baso=0, A_apical=0, 
           Aefflux=0, APTC=0, AFil=0, Aurine=0, AST=0, AabsST=0, 
           ASI=0, AabsSI=0, Afeces=0, AL=0, APlas_free=0,
           AUCCA=0, AUCCA_free=0, AUCCR=0, AUCCPTC=0, AUCfil=0, AUCKb=0,
           AUCCL=0)


################################## PBPK output #################################################################
# Time parameters
startime<- 0
stoptime<- 96
dt <- 0.1
dtout <- 0.1
Times = seq (startime, stoptime, dt)
oti = seq (1, (stoptime - startime)/ dt + 1, by = dtout/dt) # output time interval
outtime <- Times [oti]


# import data
data <- data.frame(Time = c(0, 1, 9, 15, 22, 31, 41, 50, 57, 66, 76, 85),
                   Conc = c(0, 36.8, 20.2, 18.6, 10.7, 11.2, 8.1, 7.7, 7.2, 6.08, 7.76, 4.24))



# Define the function "Modelcost" using FME package

Mcost <- function(pars) {
  
  pars = exp(pars) # Reverse to arithmetic domain from logarithmic domain
  out <- ode (y = State, times = Times, 
              parms = pars, func = pbpkmodel, method = 'lsoda') # solve the ODE equations using deSolve package
  
  outdf <- as.data.frame(out) %>% select(Time=time, Conc = CA) # select the time and serum PFOS concentrations
  
  cost <- modCost (model = outdf, obs = data, x='Time')
  
  return(cost)
}


# Select the parameter
theta.select <- log(c(
                    Free = 0.022, 
                    PL = 3.72, 
                    KbileC = 0.25
                    #Vmax_apical_invitro = 9300 
                  ))


## Markov chain Monte Carlo (MCMC) simulation ###################################
# define the configuration of parallel simulations
detectCores()                                ## check the cores
cl<- makeCluster(detectCores())              ## use all cores in our system     
registerDoParallel(cl)                       ## registers a cluster of all the cores on our system

# start time
strt<-Sys.time()

# define the iterations of MCMC simulation
niter         = 500
burninlength  = niter*0.5 ## number of initial iterations to be removed from output
outputlength  = niter*0.5 ## number of output iterations

# parallel
system.time(
  MCMC <- foreach( i = 1:4, .packages = c('magrittr','FME','truncnorm','EnvStats','invgamma','dplyr')) %dopar% {
    modMCMC(f             = Mcost,           ## defined cost function
            p             = theta.select,    ## parameters
            niter         = niter,           ## iteration number 
            jump          = 0.1,             ## jump function generate new parameters distribution with covariate matrix
            updatecov     = 5,               ## used to update the MCMC jumps.
            wvar0         = 0.1,             ## "Weight" for the initial model variance
            burninlength  = burninlength,    ## number of initial iterations to be removed from output.
            outputlength  = outputlength)    ## number of output iterations           
  }
)

#end time
print(Sys.time()-strt)

stopCluster(cl)   

## Performance four chains to check the convergences
MC.1 = as.mcmc (MCMC[[1]]$pars)   # first  chain
MC.2 = as.mcmc (MCMC[[2]]$pars)   # second chain
MC.3 = as.mcmc (MCMC[[3]]$pars)   # third  chain
MC.4 = as.mcmc (MCMC[[4]]$pars)   # fourth chain
combinedchains = mcmc.list(MC.1,MC.2,MC.3,MC.4) ## combine all chains
gelman.diag (combinedchains)          # Gelman convergence diagnosis

## Save the posterior parameters (95% CI)
quan.mouse = exp(summary(MC.1)$quantiles)  

## Trace plot using bayesplot
## Convergences plot
mcmc_trace (
  combinedchains,
  pars       = names(theta.select),
  size       = 0.5,
  facet_args = list(nrow = 2)) +
  ggplot2::scale_color_brewer()


# Best Time vs. conc. plot
out_mcmc <- ode (y = State, times = Times, parms = exp(MCMC[[1]]$bestpar), func = pbpkmodel, method = 'lsoda')
pbpkout_mcmc <- as.data.frame (out_mcmc)%>%select(Time = time, Conc=CA)

ggplot (data = pbpkout_mcmc, aes(Time, Conc)) + 
  geom_line() + geom_point(data=data, aes(Time, Conc)) 



# Create the matrix 
MC.CA   = matrix(nrow = length(Times), ncol = outputlength)

for (i in 1:outputlength) {
  
  j = i
  pars = MCMC[[1]]$pars[j,]
  out_mcmc <- ode (y = State, times = Times, 
                   parms = exp(pars), 
                   func = pbpkmodel, 
                   method = 'lsoda')
  
  MC.CA[,i] <- as.data.frame (out_mcmc)$CA
  
 
  cat("iteration = ", i , "\n") # Shows the progress of iterations, so you can see the number of iterations that has been completed, and how many left.
}

## 
MC.CA.plot <- cbind(
  Time = Times, 
  as.data.frame(t(apply(MC.CA, 1, function(y_est) c( # 1 indicates row; 
    median_est         = median(y_est,na.rm = T), 
    ci_q1              = quantile(y_est, probs = 0.25, names = FALSE,na.rm = T), 
    ci_q3              = quantile(y_est, probs = 0.75, names = FALSE,na.rm = T),
    ci_10              = quantile(y_est, probs = 0.10, names = FALSE,na.rm = T), 
    ci_90              = quantile(y_est, probs = 0.90, names = FALSE,na.rm = T),
    ci_lower_est       = quantile(y_est, probs = 0.025, names = FALSE,na.rm = T),  
    ci_upper_est       = quantile(y_est, probs = 0.975, names = FALSE,na.rm = T)  
  ))))
)

p <- 
  ggplot() + 
  geom_ribbon(data = MC.CA.plot, aes(x = Time, ymin = ci_lower_est, ymax = ci_upper_est,), 
              fill="gray") +
  geom_line(data= MC.CA.plot, aes(x = Time, y = median_est), 
            linewidth = rel(0.3), colour = "black") +
  geom_point(data = data, aes(x=Time, y= Conc), shape = 1) + theme_bw()

p

## Normalized sensitivity analysis
## Define the sensitivity function

NSC_func <- function (par) {
  
  n <- length(par) # to obtain the length of input parameters 
  m <- 7 # the model is able to output 24-hours area under curve (AUC) in seven different organs/tissues
    
  # create a n x m matrix  
  NSC  <- matrix(NA, nrow = n , ncol = 7) # this matrix contain n rows and m columns
  
  for (i in 1:n) { # create a loop to estimate the impact of individual parameter on output 
    
    New.pars   <- par %>% replace(i, as.numeric(par[i])*1.01) # replace one parameter by adding 1% of original value to the parameter 
    Rnew       <- ode (y = State, times = Times, parms = New.pars, func = pbpkmodel, method = 'lsoda') # predictions based on a 1% of change on specific parameter and other parameters keep the same parameter
    R          <- ode (y = State, times = Times, parms = par, func = pbpkmodel, method = 'lsoda') # predictions based on original parameter sets
    delta.pars <- as.numeric(par[i])/as.numeric(par[i]*0.01) # the changes from old parameter to new parameter
    
    # The estimated AUC based on new parameter set
    AUC.new  <- data.frame(Rnew) %>% filter(time == 24)%>% select(AUCCA=AUCCA,AUCCA_free=AUCCA_free, 
                                                                  AUCCR=AUCCR, AUCCPTC=AUCCPTC, AUCfil=AUCfil,
                                                                  AUCKb=AUCKb,AUCCL=AUCCL) %>% 
                                                                  select(-contains("time")) 
    # The estimated AUC based on original parameter set
    AUC.ori  <- data.frame(R)    %>% filter(time == 24)%>% select(AUCCA=AUCCA,AUCCA_free=AUCCA_free, 
                                                                  AUCCR=AUCCR, AUCCPTC=AUCCPTC, AUCfil=AUCfil,
                                                                  AUCKb=AUCKb,AUCCL=AUCCL) %>% 
                                                                  select(-contains("time"))
  
    
    ## Estimated the AUC
    for (j in 1:dim(AUC.new)[2]) {
      
      delta.AUC    =  AUC.new [,j] - AUC.ori[,j] # the AUC change from the predictions based original to new parameters

      NSC [i, j]   <- (delta.AUC/AUC.ori[,j]) * delta.pars # save the normalized sensitivity coefficients results

    }
  }
  
  rownames(NSC) = c(names(par))
  
  colnames(NSC)= c("AUCCA", "AUCCA_free", "AUCCR", 
                   "AUCCPTC", "AUCfil","AUCKb","AUCCL")
  
  return (list(NSC = NSC))
}

# Select the parameter
select.pars <- c(
  BW  = 0.4,             
  MKC = 0.0084,                                                   
  MLC = 0.034,           
  QCC = 14,             
  QLC = 0.183,           
  QKC = 0.141,          
  Htc = 0.46,            
  VPlasC =0.0312,        
  VLC = 0.035,           
  VKC = 0.0084,          
  VfilC = 8.4e-4,        
  VPTCC = 1.35e-4,       
  Free = 0.02,                      
  Vmax_baso_invitro = 393.45,        
  Km_baso = 27.2,                    
  Vmax_apical_invitro = 9300,        
  Km_apical= 52.3,                   
  RAFbaso = 4.07,                    
  RAFapi = 35,                       
  protein = 2.0e-6,                  
  GFRC = 62.1,                       
  PL = 3.72,                         
  PK = 0.8,                          
  PR = 0.2,                          
  Kdif = 0.001,                      
  Kabsc = 2.12,                                              
  KunabsC = 7.06e-5,                 
  GEC = 0.54,                        
  K0C = 1,                           
  KeffluxC = 2.49,                   
  KbileC = 0.25,                     
  KurineC = 0.32,                    
  Free = 0.022, 
  PL = 3.72, 
  KbileC = 0.25,
  Vmax_apical_invitro = 9300 
)

NSC <- NSC_func(par = select.pars)

## Heatmap
## Set the color
c2 <- colorRampPalette(c("navy", "white", "firebrick3"))(500)

# filter the zero value out 
NSC_new <- filter_if(as.data.frame(NSC), is.numeric, all_vars((.) != 0))

# make the heatmap
x<- pheatmap(as.matrix(NSC_new), cutree_rows = 3, scale = 'none',
             clustering_distance_rows = "correlation",
             cluster_cols = F,
             fontsize = 12,
             color = c2, 
             legend = T)


x


