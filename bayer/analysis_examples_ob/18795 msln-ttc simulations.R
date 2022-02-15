library(crmPack)

emptydata <- Data(doseGrid = c(1.5, 2.5, 3.5, 4.5, 6, 7))

MSLN_TTC_Model <- LogisticKadaneBayer(theta = 0.3,
                                      dmin = 1.5,
                                      dmax = 7,
                                      alpha = 1,
                                      beta = 19,
                                      shape = 0.5625,
                                      rate = 0.125)

options <- McmcOptions(burnin = 10000, step = 2, samples = 10000)
set.seed(94)

# Prior exploration
PriorSamples <- mcmc(emptydata, MSLN_TTC_Model, options)

plot(PriorSamples, MSLN_TTC_Model, emptydata) + ggtitle("Prior")

#Sample data to test Stopping Rule 1 : MTD precisely estimated: CV(MTD) <= 30%
data1 <- Data(x=c(1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 3.5, 3.5, 3.5, 4.5, 4.5, 4.5, 6, 6, 6, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5),
              y=c(  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, 0, 1, 1,   0,   0,   0,   1,   0,   0),
              cohort=c(  1,   1,   1,   2,   2,   2,   3,   3,   3,   4,   4,   4, 5, 5, 5,   6,   6,   6,   7,   7,   7),
              doseGrid=c(1.5, 2.5, 3.5, 4.5, 6, 7),
              ID=1:21)

Samples.18795.data1 <- mcmc(data1, MSLN_TTC_Model, options)
plot(Samples.18795.data1, MSLN_TTC_Model, data1)

#### ESCALATION RULES

### 1. Increments: For specifying maximum allowable increments between doses

## Rule1:
# the number of maximum dose levels to increment for the next dose
# 1, means that no dose skipping is allowed - the next dose can be maximum one level higher than the current dose.

## Rule2:
# If at least 2 DLTs out of 3 (91% probability that the dose is toxic),
# or at least 3 DLTs out of 6 (84%),
# or at least 4 DLTs out of 9 (81%) at any dose level,
# d must be lower than that dose.
# p.s. still needs to be implemented here or in NextBest class, see nextbest dose notes below relevant section.
# this is straight forward and can be calculated for any number of patients in a cohort with a user defined probability threshold of 90%
# i.e. 1- pbeta(target, x+a, n-x+b) > 0.90
# where x=DLTs, n=patients for a dose and use Beta(a=1,b=1) (or probably a=b=0.5)

#A new rule for increments must be used as crmPack only allows increase from last given dose,
#but we want to allow the increase from the maximum applied dose
Increments1 <- IncrementsNumDoseLevels(maxLevels = 1, basisLevel="max")
Increments2 <- IncrementsSafetyStopFix(toxrule = matrix(c(3,2,6,3,9,4),nrow=2))

Increments.18795 <- IncrementMin(IncrementsList=list(Increments1,Increments2))

#test Increments with data
(nextMaxDose.18795.data1   <- maxDose(Increments.18795, data1))

#NextBest.18795 <- NextBestKadaneBay(target = 0.3)
# for generation of the allocation criteria, method 2 is used in the protocol
NextBest.18795 <- NextBestMTDprob(target = 0.3, method=2)

#Test next best dose
(NextDose.18795.data1   <- nextBest(NextBest.18795, nextMaxDose.18795.data1, Samples.18795.data1, MSLN_TTC_Model, data1)$value)
(nextBest(NextBest.18795, nextMaxDose.18795.data1, Samples.18795.data1, MSLN_TTC_Model, data1)$allocation)

#Stopping rules
myStoppingfirst <- StoppingHardFirstFix(matrix(c(3,2,6,3,9,4),nrow=2), label="tox hard stop")
myStoppinglow   <- StoppingMTDdistributionBay(target = 0.3, prob = 0.8, doseeval=emptydata@doseGrid[1],
                                              direction='below', dosetested = T, label="tox first model")
myStoppinghigh  <- StoppingMTDdistributionBay(target = 0.3, prob = 0.8, doseeval=emptydata@doseGrid[emptydata@nGrid],
                                              direction='above', dosetested = T, label="safe last model")
myStoppingCV    <- StoppingMTDCV(target=0.3, thresh_cv=30)
myStoppingnpat  <- StoppingPatientsNearDose2(nPatients = 9, percentage = 0, label="already n=9")


#Stopping.18795 <- myStoppingfirst| myStoppinglow | myStoppinghigh | myStoppingCV | myStoppingnpat
Stopping.18795 <- myStoppingfirst| myStoppinglow | myStoppinghigh | myStoppingnpat

#test stopping rules with data
stopTrial(Stopping.18795, NextDose.18795.data1, Samples.18795.data1, MSLN_TTC_Model, data1)

mySize <- CohortSizeConst(3)

MSLN_TTC_Design <- BayDesign(model=MSLN_TTC_Model,
                             nextBest=NextBest.18795,
                             stopping=Stopping.18795,
                             increments=Increments.18795,
                             cohortSize=mySize,
                             data=emptydata,
                             startingDose=1.5)

set.seed(9)

examine(MSLN_TTC_Design)

# # Scenario A – Maximum dose safe
# A_Safe <- function(dose)
# {
#   MSLN_TTC_Model@prob(dose, p0=0.05, MTD=15)
# }
# # plot it in the range of the dose grid
# curve(A_Safe(x), from=1, to=8, ylim=c(0, 1))

s.a <- cbind(emptydata@doseGrid, c(0.06, 0.07, 0.08, 0.09, 0.11, 0.12))
s.b <- cbind(emptydata@doseGrid, c(0.10, 0.14, 0.21, 0.30, 0.46, 0.58))
s.c <- cbind(emptydata@doseGrid, c(0.16, 0.30, 0.50, 0.70, 0.89, 0.95))
s.d <- cbind(emptydata@doseGrid, c(0.55, 0.91, 0.99, 1.00, 1.00, 1.00))
s.e <- cbind(emptydata@doseGrid, c(0.05, 0.05, 0.05, 0.80, 0.80, 0.80))

safe  <- function(dose,scenario=s.a) {scenario[match(dose, scenario[, 1]), 2]}
late  <- function(dose,scenario=s.b) {scenario[match(dose, scenario[, 1]), 2]}
early <- function(dose,scenario=s.c) {scenario[match(dose, scenario[, 1]), 2]}
toxic <- function(dose,scenario=s.d) {scenario[match(dose, scenario[, 1]), 2]}
peak  <- function(dose,scenario=s.e) {scenario[match(dose, scenario[, 1]), 2]}

pltdat <- data.frame(dose=emptydata@doseGrid,
                     safe=safe(emptydata@doseGrid),
                     late=late(emptydata@doseGrid),
                     early=early(emptydata@doseGrid),
                     toxic=toxic(emptydata@doseGrid),
                     peak=peak(emptydata@doseGrid))

#ggplot(pltdat, aes(dose,toxic)) + geom_point() + geom_line()
ggplot(pltdat, aes(x=dose)) +
  geom_line(aes(y = safe),  color = "red") +
  geom_line(aes(y = late),  color = "steelblue") +
  geom_line(aes(y = early), color = "yellow") +
  geom_line(aes(y = toxic), color = "black") +
  geom_line(aes(y = peak),  color = "green")

nsim=2
#use.seed=819
use.seed=9
para=F

time <- system.time({safe.18795 <- simulate(MSLN_TTC_Design,
                                            args=NULL,
                                            truth=safe,
                                            nsim=nsim,
                                            seed=use.seed,
                                            mcmcOptions=options,
                                            parallel=para)

                    late.18795 <- simulate(MSLN_TTC_Design,
                                             args=NULL,
                                             truth=late,
                                             nsim=nsim,
                                             seed=use.seed,
                                             mcmcOptions=options,
                                             parallel=para)

                    early.18795 <- simulate(MSLN_TTC_Design,
                                          args=NULL,
                                          truth=early,
                                          nsim=nsim,
                                          seed=use.seed,
                                          mcmcOptions=options,
                                          parallel=para)

                    toxic.18795 <- simulate(MSLN_TTC_Design,
                                           args=NULL,
                                           truth=toxic,
                                           nsim=nsim,
                                           seed=use.seed,
                                           mcmcOptions=options,
                                           parallel=para)

                    peak.18795 <- simulate(MSLN_TTC_Design,
                                           args=NULL,
                                           truth=peak,
                                           nsim=nsim,
                                           seed=use.seed,
                                           mcmcOptions=options,
                                           parallel=para)
                    })



time[3]/60

#save.image("~/results/18795_msln_results_m3.RData")

summary(toxic.18795, toxic)
