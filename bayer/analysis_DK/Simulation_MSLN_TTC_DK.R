# You can install the development version of crmPack from github with:
# install.packages("devtools")
# devtools::install_github("Roche/crmPack")

#install and subsequently load package in R
# You can install the stable release version of crmPack from CRAN with:
# install.package("crmPack")
library(crmPack)

# install.packages("ggmcmc") to be used for plotting results produced by crmPack
library(ggmcmc)

###### NOTE ####################################################################################################
## open and run the file helpers.R in the in the crmPack installation R folder, before proceeding further below
## this is file contains Some helper functions used during development of new S4 class objects further below
## then you  need to run the newly created S$ objects / programs used below before proceeding further below
###### NOTE ####################################################################################################

## MSLN_TTC_Design Dose Escalation Simulation Set-up


## 1. model set-up
MSLN_TTC_Model <- LogisticKadane2(theta = 0.3, 
                                  dmin = 1.5, 
                                  dmax = 7, 
                                  alpha = 1, 
                                  beta = 19, 
                                  shape = 0.5625, 
                                  rate = 0.125)

## 2. escalation rules set-up
## 2.1. nextBest rule set-up
NextBestDose_Rule <- NextBestMTDCRM(target=0.3,pbomethod="max")

## 2.2. increments rules set-up 
## 2.2.1 no dose skipping 
Increment_Rule1 <-IncrementsNumDoseLevelsMaxTested(maxLevels = 1)

## 2.2.2 hard safety rule
Increment_Rule2<- IncrementsHSRBeta(target=0.3,prob=0.95,a=1,b=1)

## combine all the increments rules
Increment_Rules_Combined <- IncrementMin(IncrementsList=list(Increment_Rule1,Increment_Rule2))

## 2.3. stopping rules set-up
## 2.3.1 Minimum dose tested (1.5MBq) is toxic: Hard safety rule based on Beta distribution
# StoppingMinDoseToxicHSR <-
#   StoppingLowestDose() &
#   StoppingDoseHSRBeta(target=0.3,
#                       prob=0.95,
#                       a=1,
#                       b=1) 
# Stopping_Rule1 <- StoppingMinDoseToxicHSR



Stopping_Rule1 <- StoppingLowestDoseHSRBeta(target=0.3,
                                            prob=0.95,
                                            a=1,
                                            b=1) 

## 2.3.2 Minimum dose tested (1.5MBq) is toxic: T1 < 20%
# StoppingMinDoseToxic <-
#   StoppingLowestDose() &
#   StoppingPatientsNearDose(nPatients=3, percentage=0) &
#   StoppingTargetProb(target=c(0.3, 1), prob=0.8)
# 
# Stopping_Rule2 <- StoppingMinDoseToxic

Stopping_Rule2 <- StoppingTargetProbPatientsNearLowestDose(target=c(0.3, 1), prob=0.8,nPatients=3, percentage=0)

## 2.3.3	Maximum possible dose (7.0MBq) is safe: TN > 80% and 7.0MBq administered.
# StoppingMaxDoseSafe <-
#   StoppingHighestDose() &
#   StoppingPatientsNearDose(nPatients=3,percentage=0) &
#   StoppingTargetProb(target=c(0, 0.3),prob=0.8)

# Stopping_Rule3 <- StoppingMaxDoseSafe

Stopping_Rule3 <- StoppingTargetProbPatientsNearHighestDose(target=c(0, 0.3),prob=0.8, nPatients=3,percentage=0)

## 2.3.4 MTD precisely estimated: CV(MTD)<= 30%.
Stopping_Rule4 <- StoppingMTDCV(target = 0.30, thresh = 0.3)

## 2.3.5 Total sample size for next dose is already >= 9.
Stopping_Rule5 <- StoppingPatientsNearDose(nPatients = 9,percentage = 0)

## combine all the stopping rules  
Stopping_Rules_Combined <- (Stopping_Rule1 | Stopping_Rule2 | Stopping_Rule3  | Stopping_Rule4  | Stopping_Rule5)

## 2.3. cohortSize rules set-up
## 2.3.1 constant cohort size of 3 along the study
CohortSize_Rule1<-CohortSizeConst(3)

# combine all the cohortSize rules
CohortSize_Rules_Combined <- maxSize(CohortSize_Rule1)

## 3. data set-up
emptydata <- Data(doseGrid = c(1.5, 2.5, 3.5, 4.5, 6, 7 ))


## define simulation design   
MSLN_TTC_Design <- Design(model=MSLN_TTC_Model,
                          nextBest=NextBestDose_Rule,
                          stopping=Stopping_Rules_Combined,
                          increments=Increment_Rules_Combined,
                          cohortSize=CohortSize_Rules_Combined,
                          data=emptydata,
                          startingDose=1.5)

## define simulation seed
set.seed(94)
options <- McmcOptions(burnin=1000,step=1,samples=10000)

# Before looking at the "many trials" operating characteristics, it is important to look at
# the"single trial"operating characteristics of the dose escalation design. For this, crmPack
# provides the function examine, which generates a dataframe showing the beginning of
# several hypothetical trial courses under the design. Assuming no DLTs have been seen
# until a certain dose, then the consequences of different number of DLTs being observed
# at this dose are shown. In the current example we have

examine(object=MSLN_TTC_Design,mcmcOptions=options)

## define the true scenarios using a data generative model similar to our design


# Scenario A: Maximum dose safe 
A_Safe <- function(dose)
{
  MSLN_TTC_Model@prob(dose, p0=0.05, MTD=15)
}

# Scenario B: Late MTD
B_Late <- function(dose)
{
  MSLN_TTC_Model@prob(dose, p0=0.05, MTD=4.5)
}

# Scenario C: early MTD
C_Early <- function(dose)
{
  MSLN_TTC_Model@prob(dose, p0=0.05, MTD=2.5)
}

# Scenario D: Toxic
D_Toxic <- function(dose)
{
  MSLN_TTC_Model@prob(dose, p0=0.05, MTD=1.0)
}

# Scenario E: Mid Peak
doseProbMatrixE <- cbind(c(1.5, 2.5, 3.5, 4.5, 6, 7), c(0.05, 0.05, 0.05, 0.80, 0.80, 0.80))
E_Peak <- function(dose)
{
  doseProbMatrixE[match(dose, doseProbMatrixE[, 1]), 2]
}





# plot scenarios in the range of the dose grid
# curve(A_Safe(x), from=1, to=8, ylim=c(0, 1))
# curve(B_Late(x), from=1, to=8, ylim=c(0, 1))
# curve(C_Early(x), from=1, to=8, ylim=c(0, 1))
# curve(D_Toxic(x), from=1, to=8, ylim=c(0, 1))
# curve(E_Peak(x), from=1, to=8, ylim=c(0, 1))

# time_A_Safe <- system.time(sims_A_Safe <- simulate(MSLN_TTC_Design,
#                                                    args=NULL,
#                                                    truth=A_Safe,
#                                                    #trueMTD=NULL,
#                                                    nsim=100,
#                                                    seed=94,
#                                                    mcmcOptions=options,
#                                                    parallel=FALSE
# )
# )

time_B_Late <- system.time(sims_B_Late <- simulate(MSLN_TTC_Design,
                                                   args=NULL,
                                                   truth=B_Late,
                                                   #trueMTD=NULL,
                                                   nsim=100,
                                                   seed=94,
                                                   mcmcOptions=options,
                                                   parallel=TRUE
)
)

# time_C_Early <- system.time(sims_C_Early <- simulate(MSLN_TTC_Design,
#                                                      args=NULL,
#                                                      truth=C_Early,
#                                                      #trueMTD=NULL,
#                                                      nsim=100,
#                                                      seed=94,
#                                                      mcmcOptions=options,
#                                                      parallel=FALSE
# )
# )
# 
# time_D_Toxic <- system.time(sims_D_Toxic <- simulate(MSLN_TTC_Design,
#                                                      args=NULL,
#                                                      truth=D_Toxic,
#                                                      #trueMTD=NULL,
#                                                      nsim=100,
#                                                      seed=94,
#                                                      mcmcOptions=options,
#                                                      parallel=FALSE
# )
# )
# 
# time_E_Peak <- system.time(sims_E_Peak <- simulate(MSLN_TTC_Design,
#                                                    args=NULL,
#                                                    truth=E_Peak,
#                                                    #trueMTD=3.5,
#                                                    nsim=100,
#                                                    seed=94,
#                                                    mcmcOptions=options,
#                                                    parallel=TRUE
# )
# )

# time_A_Safe
# time_B_Late
# time_C_Early
# time_D_Toxic
# time_E_Peak
# 
# 
# sims_A_Safe
# sims_B_Late
# sims_C_Early
# sims_D_Toxic
# sims_E_Peak

# simSum_A_Safe <- summary(sims_A_Safe, truth=A_Safe)
# simSum_B_Late <- summary(sims_B_Late, truth=B_Late)
# simSum_C_Early <- summary(sims_C_Early, truth=C_Early)
# simSum_D_Toxic <- summary(sims_D_Toxic, truth=D_Toxic)
# simSum_E_Peak <- summary(sims_E_Peak, truth=E_Peak)

# simSum_A_Safe
# simSum_B_Late
# simSum_C_Early
# simSum_D_Toxic
# simSum_E_Peak

# print(plot(simSum_A_Safe))
# print(plot(simSum_B_Late))
# print(plot(simSum_C_Early))
# print(plot(simSum_D_Toxic))
# print(plot(simSum_E_Peak))


