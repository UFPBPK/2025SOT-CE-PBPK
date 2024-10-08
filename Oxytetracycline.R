library(deSolve) #load and attach add-on packages.
# time parameter
starttime <- 0
stoptime <- 24
dtout <- 0.001 # resolution of output time
Times <- seq(starttime, stoptime, dtout)
###################################################################################################################
## parameters to be changed as follows:############################################################################

# Parameters for exposure scenario
PDOSEiv = 10 # (mg/kg)
PDOSEim = 0 #(mg/kg)
PDOSEoral = 0 # (mg/kg)
#IM absorption rate constats
Kim = 0.3 # 0.15 for conventional formulation# 0.3 for long-acting formulation# IM absorption rate constant(/h)
Frac = 0.5 # 0.95 for conventional formulation#0.5 for long-acting formulation#
Kdiss = 0.02 #/h
#Dosing, multiple oral gavage
tlen = 0.001 # Length of oral gavage exposure (h/day)
tinterval = 6 # varied dependent on the exposure paradigm
tdose = 1 # dose times
Ka = 0.012 # for tablets or capsules
#Ka = 0.05 # Ka = 0.05 for experimental solution/h, intestinal absorption rate constant, 


###################################################################################################################
## Fixed parameters as follows:####################################################################################
# Physiological parameters
# Blood flow rates
QCC = 12.9 # Cardiac output (L/h/kg)
QLC = 0.297 # Fraction of flow to the liver
QKC  = 0.173 # Fraction of flow to the kidneys
QFC = 0.097 # Fraction of flow to the fat
QMC = 0.217 # Fraction of flow to the muscle
# Tissue  volumes
BW = 11.3 # Body weight(kg)
VLC = 0.0329 # Fractional liver tissue
VKC = 0.0055 # Fractional kidney tissue
VFC = 0.15 # Fractional fat tissue
VMC = 0.4565 # Fractional muscle tissue
VbloodC = 0.082 # Blood volume, fractional of BW

#Fraction of tissue volumes that is blood!
FVBF = 0.02# Blood volume fraction of fat (%)
FVBM = 0.01# Blood volume fraction of muscle (%)
FVBS = 0.01# Blood volume fraction of slowly perfused tissue (%)

# Mass Transfer Parameters (Chemical-specific parameters)
# Partition coefficients(PC, tissue:plasma)
PL = 1.89 # Liver: plasma PC
PK = 4.75 # Kidney:plasma PC
PM = 0.85# Muscle:plasma PC
PF = 0.086 #Fat:plasma PC
PR = 4.75# Richly perfused tissues:plasma PC
PS = 0.85# Slowly perfused tissues:plasma PC

#Permeability constatns (L/h/kg tissue) (Permeation area cross products)
PAFC = 0.012# Fat tissue permeability constant 
PAMC = 0.225# Muscle tissue permeability constant
PASC = 0.049# Slowly perfused tissue permeability constant

#Kinetic constats
#Oral absorbtion rate constants
Kst = 2 # /h, gastric emptying rate constant
Kint = 0.2 #/h, intestinal transit rate constant.

# IV infusion rate constants
Timeiv = 0.01 # IV injection/infusion time (h)

# Urinary elimination rate constant adjusted by bodyweight
KurineC = 0.2 # L/h/kg
# Urinary elimination rate constant
Kurine = KurineC * BW #L/h


# Cardiac output and blood flows to tissues(L/h)
QC = QCC * BW # Cardiac output
QL = QLC * QC # Liver
QK = QKC * QC # Kidney
QF = QFC * QC # Fat
QM = QMC * QC # Muscle
QR = 0.626 * QC - QK - QL # Richly perfused tissues
QS = 0.374* QC - QF - QM # Slowly perfused tissues

#  Tissue volumes (L)
VL = VLC * BW # Liver
VK = VKC * BW # Kidney
VF = VFC * BW # Fat
VM = VMC * BW # Muscle
Vblood = VbloodC * BW# Blood
VR = 0.142 * BW - VL - VK # Richly perfused tissues
VS = 0.776 * BW - VF -VM # Slowly perfused tissues

#Permeability surface area coefficients
PAF = PAFC * VF# Fat:blood permeability (L/h)
PAM = PAMC * VM# Muscle : blood permeability (L/h)
PAS = PASC * VS# Slowly perfused tissue:blood permeability (L/h)

#Volume of tissue vs blood
#Fat
VFB = FVBF * VF # fat compartment blood volume
VFT = VF - VFB# Fat compartment tissue volume

#Muscle
VMB = FVBM * VS# Muscle compartment blood volume
VMT = VM - VMB# Muslce compartment blood volume
#Slowly perfused tissue
VSB = FVBS * VS# Slowly perfused compartment blood volume
VST = VS - VSB#  Slowly perfused compartment blood volume

# Dosing
DOSEoral = PDOSEoral * BW # (mg)
DOSEiv = PDOSEiv * BW # (mg)
DOSEim = PDOSEim * BW # (mg)
# Dosing, intramuscular, dissolution model
DOSEimfast = DOSEim * Frac
DOSEimslow = DOSEim * (1 - Frac)

# IVR & oral rate constant
IVR = DOSEiv/Timeiv
oralR = DOSEoral/tlen
## PBPK model
pbpkmodel <- function(Time, State, Parmeters) {
  with(as.list(c(State, Paras)), {
    ###################################################################################################################
    ## Concentration of the chemical in vein compartment
    CVL = AL/(VL * PL) # con' of chem in liver  / PC of plasma: liver
    CVK = AK/(VK * PK) # con' of chem in Kidney / PC of plasma: kidney
    CVF = AFB/VFB
    CVR = AR/(VR * PR) # con' of chem in RPT    / PC of plasma: richly perfused tissue
    CVS = ASB/VSB # con' of chem in SPT: slowly perfused tissue
    CVM = AMB/VMB
    CMT = AMT/VMT
    CFT = AFT/VFT
    CST = ASLT/VST
    # OTCiv injection to the venous
    RIV = IVR * (Time < Timeiv)
    dAIV = RIV
    # OTC oral to the stomach
    RDOSEoral <- oralR * (Time <= tdose * tinterval) * (Time %% tinterval < tlen) 
    dADOSEoral = RDOSEoral
    RAST = RDOSEoral - Kst * AST
    dAST = RAST
    # OTCim injetion to the muscle
    RDOSEimremain = -Kdiss * DOSEimremain
    dDOSEimremain = RDOSEimremain
    Rim = Kim * Amtsite
    dAbsorb = Rim
    Rsite = - Rim + Kdiss * DOSEimremain
    dAmtsite = Rsite
    ###################################################################################################################      
    ## OTC in Intestine and Colon
    RAI = Kst * AST - Kint * AI - Ka * AI
    dAI = RAI
    Rcolon = Kint * AI
    dAcolon = Rcolon
    RAO = Ka * AI
    dAAO = RAO
    ###################################################################################################################      
    ## OTC in blood compartment
    #con' of chemical in the vein
    CV = ((QL * CVL + QK  * CVK + QF * CVF + QM * CVM + QR * CVR + QS * CVS + RIV + Rim) / QC) 
    
    CA = AA/Vblood # con' of chem in artery = amount of chem in artery / volume of blood
    RA = QC * (CV - CA) # rate of change in amount of chem in tissue of blood compartment
    dAA = RA
    dAUCCV = CV
    ###################################################################################################################    
    ## OTC in liver compartment
    RL = QL * (CA - CVL) + RAO # rate of change in amount of the chem in liver
    dAL = RL # amount of chemical in liver
    CL = AL/VL # con' of chem in liver
    dAUCCL = CL
    ###################################################################################################################    
    ## OTC in kidney compartment
    # Urinary excretion of OTC
    Rurine = Kurine * CVK
    dAurine = Rurine
    # kidney
    RK = QK * (CA - CVK) - Rurine
    dAK = RK
    CK = AK/VK # con' of chem in kieney
    dAUCCK = CK
    ###################################################################################################################
    ## OTC in muslce compartment
    RMB = QM * (CA - CVM)- PAM * CVM + (PAM * CMT)/PM
    dAMB = RMB # amount of chemical in muscle
    RMT = PAM * CVM - (PAM * CMT)/PM
    dAMT = RMT
    AMtotal = AMT + AMB
    CM =AMtotal/VM
    dAUCCM = CM
    ###################################################################################################################
    ## OTC in fat compartment
    RFT = PAF * CVF - (PAF * CFT)/PF
    dAFT = RFT
    RFB = QF * (CA - CVF) - PAF * CVF + (PAF * CFT)/PF
    dAFB = RFB # amount of chemical in fat
    AFtotal = AFT + AFB
    CF =AFtotal/VF
    ###################################################################################################################
    ## OTC in RPT of body compartment
    RR = QR * (CA - CVR)
    dAR = RR
    CR = AR/VR
    ###################################################################################################################
    ## OTC in SPT of body compartment
    RSB = QS * (CA - CVS) - PAS * CVS + (PAS * CST)/PS
    dASB = RSB
    RSLT = PAS * CVS - (PAS * CST)/PS
    dASLT = RSLT
    AStotal = ASLT + ASB
    CS =AStotal/VS
    ###################################################################################################################
    ## Mass balance
    Qbal = QC - QL - QK - QM - QF - QR - QS 
    Tmass = AA + AL + AK + Aurine + AMtotal + AFtotal + AR + AStotal
    list(c(dAbsorb, dAmtsite, dDOSEimremain, dAIV, dADOSEoral, dAST, 
           dAI, dAcolon, dAAO, dAA, dAL, dAK, dAurine, dAMB, 
           dAMT, dAUCCV, dAUCCL, dAUCCK, dAUCCM, dAFB, dAFT, 
           dAR, dASB, dASLT),c(AUCCV = AUCCV, AUCCL = AUCCL, AUCCK=AUCCK, AUCCM=AUCCM))
  })
}
State <- c(Absorb = 0, Amtsite = DOSEimfast, DOSEimremain = DOSEimslow, AIV = 0, ADOSEoral = 0, AST = 0, 
           AI = 0, Acolon = 0, AAO = 0, AA = 0, AL = 0, AK = 0,  Aurine = 0, AMB = 0, 
           AMT = 0, AUCCV = 0, AUCCL = 0, AUCCK = 0, AUCCM = 0, AFB = 0, AFT = 0, 
           AR = 0, ASB = 0, ASLT = 0)
Paras <- c(Kst, Ka, Kint, Kim, Kdiss, QC, QL, QK, QM, QF, QR, QS, PAF, PAM, PAS)

# PBPK output
out <- ode(y = State,
           times = Times, 
           func = pbpkmodel, 
           parms = Paras) 
pbpkout <- as.data.frame(out[seq(from = 1, to = dim(out)[1], by = 0.1/dtout),])
newtime <- Times[seq(from = 1, to = dim(out)[1], by = 0.1/dtout)]

## Check balance
Tmass = pbpkout$AA + pbpkout$AL + pbpkout$AK + pbpkout$Aurine + pbpkout$AMT + pbpkout$AMB + 
  pbpkout$AFT + pbpkout$AFB + pbpkout$AR + pbpkout$ASLT + pbpkout$ASB
Bal = pbpkout$AAO + pbpkout$AIV + pbpkout$Absorb - Tmass #permeability-limited model mass balance
plot(newtime, Bal, col = 'red')

## Get the data from AUC
AUCCV_dat = as.data.frame(pbpkout$AUCCV)
AUCCL_dat = as.data.frame(pbpkout$AUCCL)
AUCCK_dat = as.data.frame(pbpkout$AUCCK)
AUCCM_dat = as.data.frame(pbpkout$AUCCM)

###################################################################################################################
###################################################################################################################
###################################################################################################################
# You should change the path bellow according to your computer, and change the name in quotes of write.csv function 
# so that for each PBPK model, the name would not be replaced by latter generated files. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#setwd('Z:/AP 873/module7/files') # path 
amoumt.drug <- cbind(pbpkout, Tmass, Bal)
write.csv(amoumt.drug, 'amount of drug.csv') # change file name here for different dosing
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###################################################################################################################
###################################################################################################################
###################################################################################################################


## plot concentration in all compartments
# OTC in blood compartment
## Concentration of the chemical in vein compartment
CVL = pbpkout$AL/(VL * PL) # liver
CVK = pbpkout$AK/(VK * PK) # kidney
CVM = pbpkout$AMB/VMB # muscle
CVF = pbpkout$AFB/VFB # Fat
CVR = pbpkout$AR/(VR * PR) # richly perfused tissue
CVS = pbpkout$ASB/VSB # slowly perfused tissue
# concentration of chemical in the vein
RIV = IVR * (newtime < Timeiv)
Rim = Kim * pbpkout$Amtsite
# OTC in vein compartment
CV = ((QL * CVL + QK  * CVK + QM * CVM + QF * CVF + QR * CVR + QS * CVS + RIV + Rim) / QC) 
plot(newtime, CV, type = 'l', xlab = c('Time'),  col = 1)
# OTC in artery compartment
CA = pbpkout$AA/Vblood
plot(newtime, CA, type = 'l', xlab = c('Time'),  col = 2)

# OTC in liver compartment
CL = pbpkout$AL/VL
plot(newtime, CL, type = 'l', xlab = c('Time'),  col = 3)

# OTC in kidney compartment
CK = pbpkout$AK/VK
plot(newtime, CK, type = 'l', xlab = c('Time'),  col = 4)

# OTC in muslce compartment
AMtotal = pbpkout$AMT + pbpkout$AMB
CM = AMtotal/VM
plot(newtime, CM, type = 'l', xlab = c('Time'),  col = 5)

# OTC in fat compartment
AFtotal = pbpkout$AFT + pbpkout$AFB
CF = AFtotal/VF
plot(newtime, CF, type = 'l', xlab = c('Time'),  col = 6)

# OTC in RPT compartment
CR = pbpkout$AR/VR
plot(newtime, CR, type = 'l', xlab = c('Time'),  col = 7)

# OTC in SPT compartment
AStotal = pbpkout$ASLT + pbpkout$ASB
CS = AStotal/VS
plot(newtime, CS, type = 'l', xlab = c('Time'),  col = 8)

###################################################################################################################
###################################################################################################################
###################################################################################################################
# you should change the path bellow according to your computer, and change the name in quotes of write.csv function 
# so that for each PBPK model, the name would not be replaced by latter generated files. 
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
concentration.durg <- data.frame(newtime, CV, CA, CL, CK, CM, CF, CR, CS, CVL, CVK, CVM, CVF, CVR, CVS)
write.csv(concentration.durg, 'concentration of drug-XXmg.csv') # change file name here for different dosing
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###################################################################################################################
###################################################################################################################
###################################################################################################################


