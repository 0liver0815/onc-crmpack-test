#Code used in:
#Journal of Statistical Software
#May 2019, Volume 89, Issue 10
#Model-Based Dose Escalation Designs in R with crmPack

library(crmPack)
coarseGrid <- c(25, 50, 100, 200, 300)
model <- MinimalInformative(dosegrid = coarseGrid, refDose = 100,
                            logNormal = TRUE, threshmin = 0.1, threshmax = 0.2, seed = 432,
                            control = list(max.time = 30))$model
PL <- 0.001
data <- Data(x = c(PL, 25, 25, 25, PL, 50, 50, 50, PL, 100, 100, 100),
             y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
             cohort = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
             doseGrid = c(PL, seq(25, 300, 25)), ID = 1:12, placebo = TRUE)

plot(data)
plot(data, blind = TRUE)

options <- McmcOptions(burnin = 1000, step = 2, samples = 10000)
set.seed(94)
samples <- mcmc(data, model, options)

plot(samples, model, data) + ggtitle("Posterior")
emptydata <- Data(doseGrid = data@doseGrid, placebo = TRUE)
priorsamples <- mcmc(emptydata, model, options)
plot(priorsamples, model, emptydata) + ggtitle("Prior")

myIncrements <- IncrementsRelative(intervals = c(0, 100, 200),
                                   increments = c(1, 0.5, 0.33))

(nextMaxDose <- maxDose(myIncrements, data))

myNextBest <- NextBestNCRM(target = c(0.2, 0.35), overdose = c(0.35, 1),
                           maxOverdoseProb = 0.25)

nextDoseRes <- nextBest(myNextBest, nextMaxDose, samples, model, data)
(nextDoseVal <- nextDoseRes$value)

(nextDoseRes$plot)

myStopping1 <- StoppingMinPatients(nPatients = 30)
myStopping2 <- StoppingTargetProb(target = c(0.2, 0.35), prob = 0.5)
myStopping3 <- StoppingPatientsNearDose(nPatients = 9, percentage = 20)
myStopping  <- myStopping1 | (myStopping2 & myStopping3)

stopTrial(myStopping, nextDoseVal, samples, model, data)


mySize <- CohortSizeConst(3)
mySizePL <- CohortSizeConst(1)
design <- Design(model = model, nextBest = myNextBest,
                 stopping = myStopping, increments = myIncrements, cohortSize = mySize,
                 PLcohortSize = mySizePL, data = emptydata, startingDose = 25)

set.seed(23)
examine(design, options)

myTruth <- function(dose) {
  model@prob(dose, alpha0 = 4.5, alpha1 = 8)
}

doseProbMatrix <- cbind(c(1, 2, 3, 4, 5), c(0.01, 0.02, 0.04, 0.06, 0.09))
myTruthMatrix <- function(dose) {
  doseProbMatrix[match(dose, doseProbMatrix[, 1]), 2]
}

#mySimsTime <- system.time({
#  mySims <- simulate(design, truth = myTruth, nsim = 100,
#                     seed = 819, mcmcOptions = options, parallel = FALSE)})[3]

mySimsTime_p <- system.time({
  mySims <- simulate(design, truth = myTruth, nsim = 100,
                     seed = 819, mcmcOptions = options, parallel = TRUE)})[3]

# mySimsTime_1000 <- system.time({
#   mySims <- simulate(design, truth = myTruth, nsim = 1000,
#                      seed = 819, mcmcOptions = options, parallel = FALSE)})[3]

mySimsTime_p1000 <- system.time({
  mySims <- simulate(design, truth = myTruth, nsim = 1000,
                     seed = 819, mcmcOptions = options, parallel = TRUE)})[3]

