## install.packages 
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

# Physiological parameters
# Body weight and fraction of blood flow to tissues
BW    = 0.273      # kg,       body weight                       ; Measured data from Liang et al., 2021
QCC   = 15.0       # L/h/kg^0.75, Cardiac output; Value obtained ; Brown et al., 1997
QLuC  = 1          # % of QCC, Fraction of blood flow to lung    ; Brown et al., 1997
QGIC  = 0.101      # % of QCC, Fraction of blood flow to spleen  ; Davis and Morris, 1993
QLC   = 0.021      # % of QCC, Fraction of blood flow to Liver   ; Brown et al., 1997
QKC   = 0.141      # % of QCC, Fraction of blood flow to kidney  ; Brown et al., 1997
QSC   = 0.0085     # % of QCC, Fraction of blood flow to brain   ; Brown et al., 1997; Upton, 2008

# Fraction of volume of tissues out of body weight
VLuC  = 0.0074     # % of BW,  Fraction of volume of lung        ; Measured data from Dr. Kreyling
VGIC  = 0.1014     # % of BW,  Fraction of volume of GI          ; Measured data from Dr. Kreyling
VLC   = 0.0346     # % of BW,  Fraction of volume of liver       ; Measured data from Dr. Kreyling
VKC   = 0.0091     # % of BW,  Fraction of volume of kidney      ; Measured data from Dr. Kreyling
VSC   = 0.0038     # % of BW,  Fraction of volume of spleen      ; Measured data from Dr. Kreyling
VBC   = 0.0717     # % of BW,  Fraction of volume of blood       ; Measured data from Dr. Kreyling

# Fraction of blood volume in tissues
BVLu   = 0.36      # % of VLu, Fraction of blood volume in lung  ; Brown et al., 1997 
BVGI   = 0.04      # % of VGI, Fraction of blood volume in GI    ; Estimated from MacCalman and Tran, 2009
BVL    = 0.21      # % of VL,  Fraction of blood volume in liver ; Brown et al., 1997
BVK    = 0.16      # % of VLK, Fraction of blood volume in kidney; Brown et al., 1997
BVS    = 0.22      # % of VLS, Fraction of blood volume in spleen; Brown et al., 1997
BVR    = 0.04      # % of VR,  Fraction of blood volume in rest of body; Brown et al., 1997

# Partition coefficient
PLu    = 0.15      # Unitless, Partition coefficient of lung     ; Lin et al., 2016
PGI    = 0.15      # Unitless, Partition coefficient of GI       ; Lin et al., 2016
PL     = 0.08      # Unitless, Partition coefficient of Liver    ; Lin et al., 2016
PK     = 0.15      # Unitless, Partition coefficient of Kidney   ; Lin et al., 2016
PS     = 0.15      # Unitless, Partition coefficient of Spleen   ; Lin et al., 2016
PR     = 0.15      # Unitless, Partition coefficient of rest of body; Lin et al., 2016

# Membrane-limited permeability
PALuC  = 0.001     # Unitless, Membrane-limited permeability coefficient of lung         ; Lin et al., 2016
PAGIC  = 0.001     # Unitless, Membrane-limited permeability coefficient of GI           ; Lin et al., 2016
PALC   = 0.001     # Unitless, Membrane-limited permeability coefficient of Liver        ; Lin et al., 2016
PAKC   = 0.001     # Unitless, Membrane-limited permeability coefficient of Kidney       ; Lin et al., 2016 
PASC   = 0.001     # Unitless, Membrane-limited permeability coefficient of Spleen       ; Lin et al., 2016
PARC   = 0.001     # Unitless, Membrane-limited permeability coefficient of rest of body ; Lin et al., 2016

# Endocytic parameters in liver; Chou et al., 2023_Table S2
KLRESrelease      = 10       # 1/h,      Release rate constant of Phagocytic cells 
KLRESmax          = 20       # 1/h,      Maximum uptake rate constant of phagocytic cells
KLRES50           = 24       # h,        Time reaching half maximum uptake rate
KLRESn            = 0.5      # Unitless, Hill coefficient, (unitless)
ALREScap          = 100      # ug/g,     tissue, Uptake capacity per tissue weight 

# Endocytic parameters in spleen
KSRESrelease      = 1        # 1/h,      Release rate constant of Phagocytic cells
KSRESmax          = 5        # 1/h,      Maximum uptake rate constant of phagocytic cells
KSRES50           = 24       # h,        Time reaching half maximum uptake rate
KSRESn            = 0.5      # Unitless, Hill coefficient, (unitless)
ASREScap          = 100      # ug/g,     tissue, Uptake capacity per tissue weight

# Endocytic parameters in kidney
KKRESrelease      =  1       # 1/h,      Release rate constant of phagocytic cells
KKRESmax          =  5       # 1/h,      Maximum uptake rate constant of phagocytic cells
KKRES50           =  24      # h,        Time reaching half maximum uptake rate
KKRESn            =  0.5     # Unitless, Hill coefficient, (unitless)
AKREScap          =  100     # ug/g,     tissue, Uptake capacity per tissue weight

# Endocytic parameters in lung
KpulRESrelease     = 1       # 1/h,      Release rate constant of phagocytic cells
KpulRESmax         = 5       # 1/h,      Maximum uptake rate constant of phagocytic cells
KpulRES50          = 24      # h,        Time reaching half maximum uptake rate
KpulRESn           = 0.5     # Unitless, Hill coefficient, (unitless)
ApulREScap         = 100     # ug/g,     tissue, Uptake capacity per tissue weight

#+ Endocytic parameters in interstitium of lung
KinRESrelease      = 7       # 1/h,      Release rate constant of phagocytic cells
KinRESmax          = 20      # 1/h,      Maximum uptake rate constant of phagocytic cells
KinRES50           = 24      # h,        Time reaching half maximum uptake rate
KinRESn            = 0.5     # Unitless, Hill coefficient, (unitless)
AinREScap          = 100     # ug/g,     tissue, Uptake capacity per tissue weight

# Endocytic parameters in rest of body
KRRESrelease       = 0.2     # 1/h,      Release rate constant of phagocytic cells
KRRESmax           = 0.1     # 1/h,      Maximum uptake rate constant of phagocytic cells
KRRES50            = 24      # h,        Time reaching half maximum uptake rate
KRRESn             = 0.5     # Unitless, Hill coefficient, (unitless)
ARREScap           = 100     # ug/g,     tissue, Uptake capacity per tissue weight

# Endocytic parameters in GI tract
KGIRESrelease      = 0.001   # 1/h,      Release rate constant of phagocytic cells
KGIRESmax          = 0.1     # 1/h,      Maximum uptake rate constant of phagocytic cells
AGIREScap          = 100     # ug/g,     tissue, Uptake capacity per tissue weight

# Uptake and elimination parameters
KGIb               = 5.41e-3  # 1/h, absorption rate of GI tract
KbileC             = 0.00003  # L/h/kg^0.75, Biliary clearance
KfecesC            = 0.141    # L/h/kg^0.75, Fecal clearance
KurineC            = 0.00003  # L/h/kg^0.75, Urinary clearance

#+ Respiratory tract compartment--Membrane-limited model with phagocytic cells in pulmonary and interstitial regions 
#+ Transport or clearance rate in lung (1/h)
Kuab               = 8.36e-5   # 1/h, Transport rate constant from upper airway to brain
Kpulrestra         = 8.65e-4   # 1/h, Transport rate constant from inactive pulmonary phagocytizing cells to trachealbronchial region
Kpulin             = 0.126     # 1/h, Transport rate constant from pulmonary region to interstitium of lung
Kinpul             = 1.12e-6   # 1/h, Transport rate constant from interstitium of lung to pulmonary region
Ktragi             = 5.52e-3   # 1/h, Transport rate constant from trachealbronchial region to the GI lumen
Kuagi              = 0.335     # 1/h, Transport rate constant from upper airways to the GI lumen

#Parameters for exposure scenario
PDOSEoral = 0.1                # (mg)

# Blood flow to tissues (L/h)
QC     = QCC*BW^0.75
QGI    = QC*QGIC
QL     = QC*QLC
QK     = QC*QKC
QS     = QC*QSC
QR     = QC*(1 - QGIC - QLC - QKC - QSC)

# Tissue volumes (L)
VLu    = BW*VLuC
VGI    = BW*VGIC
VL     = BW*VLC
VK     = BW*VKC
VS     = BW*VSC
VB     = BW*VBC
VR     = BW*(1 - VLuC - VGIC - VLC - VKC - VSC - VBC)  
VBal   = BW - (VLu + VGI + VL + VK + VS + VR + VB)

# Tissue volumes for different compartments (L)
VLub   = VLu*BVLu  
VLut   = VLu-VLub  
VGIb   = VGI*BVGI 
VGIt   = VGI-VGIb
VLb    = VL*BVL  
VLt    = VL-VLb
VKb    = VK*BVK    
VKt    = VK-VKb    
VSb    = VS*BVS   
VSt    = VS-VSb
VRb    = VR*BVR   
VRt    = VR-VRb   

# Permeability coefficient-surface area cross-product (L/h)
PALu   = PALuC*QC  
PAGI   = PAGIC*QGI 
PAL    = PALC*QL 
PAK    = PAKC*QK 
PAS    = PASC*QS
PAR    = PARC*QR 


# Excretion 
Kbile       = KbileC*BW^0.75
Kurine      = KurineC*BW^0.75
Kfeces      = KfecesC*BW^0.75

#Dosing
IVDOSE1 = 32.9352  #(mg)


#State represents a vector of initial values of A; Parameters are K and PS, etc.
pbpkmodel <- function(Time, State, Parameters) {
  with(as.list(State), { 
      with (as.list(Parameters), {
    ##---------------------------------------------------------------------##
    ##Concentration of the chemical in the vein or tissue of each compartment
    CA        = AA/ (VB*0.2)
    CV        = AV/ (VB*0.8)
    Cin       = Ain/VLut                      
    CVLu      = ALub/VLub
    CVS       = ASb/VSb
    CSt       = ASt/VSt
    CVGI      = AGIb/VGIb
    CGIt      = AGIt/VGIt
    CVL       = ALb/VLb
    CLt       = ALt/VLt
    CVK       = AKb/VKb
    CKt       = AKt/VKt
    CVR       = ARb/VRb 
    CRt       = ARt/VRt
    
    # Endocytosis rate (1/h)
    KLRESUP     = KLRESmax*(1-(ALRES/(ALREScap*VL)))
    KSRESUP     = KSRESmax*(1-(ASRES/(ASREScap*VS)))
    KKRESUP     = KKRESmax*(1-(AKRES/(AKREScap*VK)))
    KpulRESUP   = KpulRESmax*(1-(ApulRES/(ApulREScap*VLu)))
    KinRESUP    = KinRESmax*(1-(AinRES/(AinREScap*VLu)))
    KGIRESUP    = KGIRESmax*(1-(AGIRES/(AGIREScap*VGI)))
    KRRESUP     = KRRESmax*(1-(ARRES/(ARREScap*VR)))
    
    #+ Equation for estimation of the rate of each compartment 
    RAA       = QC*CVLu - QC*CA                                             
    RAV       = QL*CVL + QK*CVK + QR*CVR - QC*CV                               
    RARb      = QR*(CA - CVR) - PAR*CVR + (PAR*CRt)/PR
    RARt      = PAR*CVR - (PAR*CRt)/PR - (KRRESUP*ARt - KRRESrelease*ARRES) + Kuab*Aua
    RARRES    = KRRESUP*ARt - KRRESrelease*ARRES
    RAua      = -(Kuagi + Kuab)*Aua
    Rtra      = Kpulrestra*ApulRES - Ktragi*Atra
    RApul     = Kinpul*Ain - Kpulin*Apul - KpulRESUP*Apul + KpulRESrelease*ApulRES
    RApulRES  = KpulRESUP*Apul - KpulRESrelease*ApulRES - Kpulrestra*ApulRES
    RAin      = Kpulin*Apul - Kinpul*Ain + PALu*CVLu - (PALu*Cin)/PLu - KinRESUP*Ain + KinRESrelease*AinRES
    RAinRES   = KinRESUP*Ain - KinRESrelease*AinRES
    RALub     = QC*(CV-CVLu) - PALu*CVLu + (PALu*Cin)/PLu
    RASb      = QS*(CA-CVS) - PAS*CVS + (PAS*CSt)/PS
    RASt      = PAS*CVS - (PAS*CSt)/PS - KSRESUP*ASt + KSRESrelease*ASRES
    RASRES    = KSRESUP*ASt-KSRESrelease*ASRES
    RAGIb     = QGI*(CA-CVGI) - PAGI*CVGI + (PAGI*CGIt)/PGI + KGIb*ALumen
    RAGIt     = PAGI*CVGI - (PAGI*CGIt)/PGI - (KGIRESUP*AGIt-KGIRESrelease*AGIRES)
    RALumen   = Kuagi*Aua + Ktragi*Atra + Kbile*CLt - (Kfeces + KGIb)*ALumen
    RAGIRES   = KGIRESUP*AGIt - KGIRESrelease*AGIRES
    RALb      = QL*(CA-CVL) + (QS*CVS) + (QGI*CVGI) - PAL*CVL + (PAL*CLt)/PL - KLRESUP*ALb + KLRESrelease*ALRES
    RALt      = PAL*CVL - (PAL*CLt)/PL - Kbile*CLt
    RALRES    = KLRESUP*ALb - KLRESrelease*ALRES
    RKb       = QK*(CA-CVK) - PAK*CVK + (PAK*CKt)/PK - Kurine*CVK
    RKt       = PAK*CVK - (PAK*CKt)/PK - KKRESUP*AKt + KKRESrelease*AKRES
    RKRES     = KKRESUP*AKt-KKRESrelease*AKRES
    Rurine    = Kurine*CVK
    Rbile     = Kbile*CLt
    Rfeces    = Kfeces*ALumen

    #+ ODE equations for compartments in the female rat
    dAA          = RAA
    dAV          = RAV
    dARb         = RARb
    dARt         = RARt
    dARRES       = RARRES
    dAua         = RAua
    dAtra        = Rtra
    dApul        = RApul
    dApulRES     = RApulRES
    dAin         = RAin
    dAinRES      = RAinRES
    dALub        = RALub
    dASb         = RASb
    dASt         = RASt
    dASRES       = RASRES
    dAGIb        = RAGIb
    dAGIt        = RAGIt
    dALumen      = RALumen
    dAGIRES      = RAGIRES
    dALb         = RALb
    dALt         = RALt
    dALRES       = RALRES
    dAKb         = RKb
    dAKt         = RKt
    dAKRES       = RKRES
    dAurine      = Rurine
    dAbile       = Rbile
    dAfeces      = Rfeces

  #+ Total amount of NPs in tissues
  ABlood    = AA + AV;
  ALung     = Atra + Apul + Ain + ALub + ApulRES + AinRES;
  ASpleen   = ASb + ASt + ASRES;
  ARest     = ARb + ARt + ARRES;
  AGI       = AGIb + AGIt + ALumen + AGIRES;
  ALiver    = ALb + ALt + ALRES;
  AKidney   = AKb + AKt + AKRES;

  #+ Total concentrations of NPs in Tissues
  CBlood    = ABlood/VB;
  CLung     = ALung/VLu;
  CSPleen   = ASpleen/VS;
  CRest     = ARest/VR;
  CGI       = AGI/VGI;
  CLiver    = ALiver/VL;
  CKidney   = AKidney/VK;

  #+ AUC
  dAUCB        = CBlood;
  dAUCLu       = CLung;
  dAUCS        = CSPleen;
  dAUCRt       = CRest;
  dAUCGI       = CGI;
  dAUCL        = CLiver;
  dAUCK        = CKidney;


    
    QBal   = QC - (QGI + QL + QK + QS + QR)
    Tmass  = ABlood + Aua + ALung + AGI + ALiver + ASpleen + AKidney + ARt + Aurine + Afeces
    BAL    = IVDOSE1 - Tmass
    
    
    list(c(dAA, dAV, dARb, dARt, dARRES, dAua, dAtra, dApul, dApulRES, dAin, dAinRES,
           dALub, dASb, dASt, dASRES, dAGIb, dAGIt, dALumen, dAGIRES, dALb, dALt,
           dALRES, dAKb, dAKt, dAKRES, dAurine, dAbile, dAfeces, dAUCB, dAUCLu, dAUCS,
           dAUCRt, dAUCGI, dAUCL, dAUCK), c(QBal = QBal))
      
    }) # end with parameters
  }) # end with y
} # end bd.model

  
  
# mg, Amount of MPs in tissues/capillary blood/phagocytic cells of compartment i
State <- c(AA = 0, AV = 0, ARb = 0, ARt = 0, ARRES = 0, Aua = 0, Atra = 0, Apul = 0, ApulRES=0,
           Ain =0, AinRES=0, ALub = 0, ASb =0, ASt=0, ASRES=0, AGIb=0, AGIt=0, ALumen=0, AGIRES=0,
           ALb=0, ALt = 0, ALRES=0, AKb=0, AKt=0, AKRES=0, Aurine=0,Abile=0, Afeces=0, AUCB=0,
           AUCLu=0, AUCS=0, AUCRt=0, AUCGI=0, AUCL=0, AUCK=0)

Parameters <- c()

#Time parameters
starttime <- 0
stoptime <- 24
dt <- 0.1 
Times <- seq(starttime, stoptime, dt)

#PBPK output, solve the differential eq. w/ parameters
out <- ode(y= State,
           times = Times,
           func = pbpkmodel,
           parms = Parameters,
           method = 'lsoda')

pbpkout <- as.data.frame(out)

## select parameters
theta.select <- c(
  # Body weight and fraction of blood flow to tissues
  BW    = 0.273,     # kg,       body weight                       ; Measured data from Liang et al., 2021
  QCC   = 15.0 ,     # L/h/kg^0.75, Cardiac output; Value obtained ; Brown et al., 1997
  QLuC  = 1,         # % of QCC, Fraction of blood flow to lung    ; Brown et al., 1997
  QGIC  = 0.101,     # % of QCC, Fraction of blood flow to spleen  ; Davis and Morris, 1993
  QLC   = 0.021,     # % of QCC, Fraction of blood flow to Liver   ; Brown et al., 1997
  QKC   = 0.141,     # % of QCC, Fraction of blood flow to kidney  ; Brown et al., 1997
  QSC   = 0.0085,    # % of QCC, Fraction of blood flow to brain   ; Brown et al., 1997; Upton, 2008
  
  
  # Partition coefficient
  PLu    = 0.15,      # Unitless, Partition coefficient of lung     ; Lin et al., 2016
  PGI    = 0.15,      # Unitless, Partition coefficient of GI       ; Lin et al., 2016
  PL     = 0.08,      # Unitless, Partition coefficient of Liver    ; Lin et al., 2016
  PK     = 0.15,      # Unitless, Partition coefficient of Kidney   ; Lin et al., 2016
  PS     = 0.15,      # Unitless, Partition coefficient of Spleen   ; Lin et al., 2016
  PR     = 0.15 ,     # Unitless, Partition coefficient of rest of body; Lin et al., 2016
  
  # Membrane-limited permeability
  PALuC  = 0.001,     # Unitless, Membrane-limited permeability coefficient of lung         ; Lin et al., 2016
  PAGIC  = 0.001 ,    # Unitless, Membrane-limited permeability coefficient of GI           ; Lin et al., 2016
  PALC   = 0.001,     # Unitless, Membrane-limited permeability coefficient of Liver        ; Lin et al., 2016
  PAKC   = 0.001,     # Unitless, Membrane-limited permeability coefficient of Kidney       ; Lin et al., 2016 
  PASC   = 0.001,     # Unitless, Membrane-limited permeability coefficient of Spleen       ; Lin et al., 2016
  PARC   = 0.001,     # Unitless, Membrane-limited permeability coefficient of rest of body ; Lin et al., 2016
  
  
  
  # Endocytic parameters in spleen
  KSRESrelease      = 1,         # 1/h,      Release rate constant of phagocytic cells
  KSRESmax          = 5,         # 1/h,      Maximum uptake rate constant of phagocytic cells
  KSRES50           = 24,        # h,        Time reaching half maximum uptake rate
  KSRESn            = 0.5,       # Unitless, Hill coefficient, (unitless)
  ASREScap          = 100,       # ug/g,     tissue, Uptake capacity per tissue weight
  
  # Endocytic parameters in kidney
  KKRESrelease      =  1,        # 1/h,      Release rate constant of phagocytic cells
  KKRESmax          =  5,        # 1/h,      Maximum uptake rate constant of phagocytic cells
  KKRES50           =  24,       # h,        Time reaching half maximum uptake rate
  KKRESn            =  0.5,      # Unitless, Hill coefficient, (unitless)
  AKREScap          =  100,      # ug/g,     tissue, Uptake capacity per tissue weight
  
  # Endocytic parameters in lung
  KpulRESrelease     = 1,        # 1/h,      Release rate constant of phagocytic cells
  KpulRESmax         = 5,        # 1/h,      Maximum uptake rate constant of phagocytic cells
  KpulRES50          = 24,       # h,        Time reaching half maximum uptake rate
  KpulRESn           = 0.5,      # Unitless, Hill coefficient, (unitless)
  ApulREScap         = 100,      # ug/g,     tissue, Uptake capacity per tissue weight
  
  #+ Endocytic parameters in interstitium of lung
  KinRESrelease      = 7,        # 1/h,      Release rate constant of phagocytic cells
  KinRESmax          = 20,       # 1/h,      Maximum uptake rate constant of phagocytic cells
  KinRES50           = 24,       # h,        Time reaching half maximum uptake rate
  KinRESn            = 0.5 ,     # Unitless, Hill coefficient, (unitless)
  AinREScap          = 100,      # ug/g,     tissue, Uptake capacity per tissue weight
  
  # Endocytic parameters in rest of body
  KRRESrelease       = 0.2,      # 1/h,      Release rate constant of phagocytic cells
  KRRESmax           = 0.1,      # 1/h,      Maximum uptake rate constant of phagocytic cells
  KRRES50            = 24,       # h,        Time reaching half maximum uptake rate
  KRRESn             = 0.5,      # Unitless, Hill coefficient, (unitless)
  ARREScap           = 100,      # ug/g,     tissue, Uptake capacity per tissue weight
  
  # Endocytic parameters in GI tract
  KGIRESrelease      = 0.001,    # 1/h,      Release rate constant of phagocytic cells
  KGIRESmax          = 0.1,      # 1/h,      Maximum uptake rate constant of phagocytic cells
  AGIREScap          = 100       # ug/g,     tissue, Uptake capacity per tissue weight
)


## Normalized sensitivity analysis
## Define the sensitivity function
NSC_func <- function (par, State) {
  
  n <- length(par)
  # create a matrix  
  NSC  <- matrix(NA, nrow = n , ncol = 6) 
  
  for (i in 1:n) {
    
    New.pars   <- par %>% replace(i, as.numeric(par[i])*1.01)
    Rnew       <- ode (y = State, times = Times, parms = New.pars, func = pbpkmodel, method = 'lsoda')
    R          <- ode (y = State, times = Times, parms = par, func = pbpkmodel, method = 'lsoda')
    delta.pars <- as.numeric(par[i])/as.numeric(par[i]*0.01)
    
    AUC.new  <- data.frame(Rnew) %>% filter(time == 24)%>% select(AUCB=AUCB, AUCLu=AUCLu,
                                                                  AUCS=AUCS,AUCRt=AUCRt,
                                                                  AUCL=AUCL,AUCK=AUCK) %>% 
                                                                  select(-contains("time")) # select the AUC i
    
    AUC.ori  <- data.frame(R)    %>% filter(time == 24)%>% select(AUCB=AUCB, AUCLu=AUCLu,
                                                                  AUCS=AUCS,AUCRt=AUCRt,
                                                                  AUCL=AUCL,AUCK=AUCK)  %>% 
                                                                  select(-contains("time"))
    
    
    ## Estimated the AUC
    for (j in 1:dim(AUC.new)[2]) {
      
      delta.AUC    =  AUC.new [,j] - AUC.ori[,j]
      
      NSC [i, j]   <- (delta.AUC/AUC.ori[,j]) * delta.pars
      
    }
  }
  
  rownames(NSC) = c(names(par))
  
  colnames(NSC)= c("AUCB", "AUCLu", "AUCS", "AUCRt", "AUCL","AUCK")
  
  return (NSC)
}

State_iv = State
State_oral = State
State_iv['AV'] = 32.935 
State_oral["ALumen"]= 32.935 

NSC.IV <- NSC_func(par = theta.select, State = State_iv)
NSC.oral <- NSC_func(par = theta.select, State = State_oral)


## Heatmap
## Set the color
c2 <- colorRampPalette(c("navy", "white", "firebrick3"))(500)


NSC_new_iv <- filter_if(as.data.frame(NSC.IV), is.numeric, all_vars((.) != 0))
NSC_new_oral <- filter_if(as.data.frame(NSC.oral), is.numeric, all_vars((.) != 0))


# IV exposure
 pheatmap(as.matrix(NSC_new_iv), cutree_rows = 3, scale = 'row',
               clustering_distance_rows = "correlation",
               cluster_cols = F,
               fontsize = 12,
               color = c2, 
               legend = T)
  
# Oral exposure  
 pheatmap(as.matrix(NSC_new_oral), cutree_rows = 3, scale = 'row',
                clustering_distance_rows = "correlation",
                cluster_cols = F,
                fontsize = 12,
                color = c2, 
                legend = T)
  
  


