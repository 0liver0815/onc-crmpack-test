#Code used in:
#Journal of Statistical Software
#May 2019, Volume 89, Issue 10
#Model-Based Dose Escalation Designs in R with crmPack

library(crmPack)

.OneParExp <- setClass(Class = "OneParExp",
                       contains = "Model",
                       representation(skeletonFun = "function",
                                      skeletonProbs = "numeric",
                                      lambda = "numeric"))

OneParExp <- function(skeletonProbs, doseGrid, lambda) {
  skeletonFun <- approxfun(x = doseGrid, y = skeletonProbs, rule = 2)
  invSkeletonFun <- approxfun(x = skeletonProbs, y = doseGrid, rule = 1)
  .OneParExp(skeletonFun = skeletonFun, skeletonProbs = skeletonProbs, lambda = lambda,
             datamodel = function() {
               for (i in 1:nObs) {
                 y[i] ~ dbern(p[i])
                 p[i] <- skeletonProbs[xLevel[i]]^theta
               }
             },
             datanames = c("nObs", "y", "xLevel"),
             prob = function(dose, theta) {
               skeletonFun(dose)^theta
             },
             dose = function(prob, theta) {
               invSkeletonFun(prob^(1 / theta))
             },
             priormodel = function() {
               theta ~ dexp(lambda)
             },
             modelspecs = function() {
               list(skeletonProbs = skeletonProbs, lambda = lambda)
             },
             init = function() {
               list(theta = 1)
             }, sample = "theta")
}

#code from section 3.1 to run example
options <- McmcOptions(burnin = 1000, step = 2, samples = 10000)
PL <- 0.001
data <- Data(x = c(PL, 25, 25, 25, PL, 50, 50, 50, PL, 100, 100, 100),
             y = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0),
             cohort = c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
             doseGrid = c(PL, seq(25, 300, 25)), ID = 1:12, placebo = TRUE)
myIncrements <- IncrementsRelative(intervals = c(0, 100, 200),
                                   increments = c(1, 0.5, 0.33))
(nextMaxDose <- maxDose(myIncrements, data))
########################

(skeletonProbs <- round(data@doseGrid / max(data@doseGrid) / 2, 2))

newModel <- OneParExp(skeletonProbs = skeletonProbs,
                      doseGrid = data@doseGrid, lambda = 1)

newSamples <- mcmc(data, newModel, options)
plot(newSamples, newModel, data)

.NextBestMinDist <- setClass(Class = "NextBestMinDist",
                             contains = "NextBest", representation(target = "numeric"))

NextBestMinDist <- function(target) {
  .NextBestMinDist(target = target)
}

setMethod("nextBest",
          signature = signature(nextBest = "NextBestMinDist",
                                doselimit = "numeric", samples = "Samples", model = "Model",
                                data = "Data"),
          def = function(nextBest, doselimit, samples, model, data, ...) {
            dosesOK <- if (length(doselimit)) {
              which(data@doseGrid <= doselimit)
            } else {
              seq_along(data@doseGrid)
            }
            modelfit <- fit(samples, model, data)
            probDLT <- modelfit$middle[dosesOK]
            doses <- modelfit$dose[dosesOK]
            bestIndex <- which.min(abs(probDLT - nextBest@target))
            bestDose <- doses[bestIndex]
            return(list(value = bestDose))
          })

newMyNextBest <- NextBestMinDist(target = 0.3)
(newNextDoseVal <- nextBest(newMyNextBest, nextMaxDose, newSamples,
                            newModel, data)$value)
