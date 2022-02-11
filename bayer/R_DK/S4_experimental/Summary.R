

################################################################################

##' Class for the summary of model-based simulations output
##'
##' In addition to the slots in the parent class
##' \code{\linkS4class{GeneralSimulationsSummary}}, it contains two slots with
##' model fit information.
##'
##' Note that objects should not be created by users, therefore no
##' initialization function is provided for this class.
##'
##' @slot fitAtDoseMostSelected fitted toxicity rate at dose most often selected
##' @slot meanFit list with the average, lower (2.5%) and upper (97.5%)
##' quantiles of the mean fitted toxicity at each dose level
##'
##' @export
##' @keywords classes
.SimulationsSummary <-
  setClass(Class="SimulationsSummary",
           representation(fitAtDoseMostSelected="numeric",
                          meanFit="list",
                          trialMTDaboveTopDose="numeric",
                          trialMTDbelowFirstDose="numeric",
                          medianMTD="numeric",
                          meanCVMTD="numeric",
                          medianMTDRE="numeric",
                          nCohorts="ANY",
                          nOD="ANY",
                          trueMTD="ANY",
                          stopReasonsUniquePercent="numeric",
                          stopReasonsFirstHitFreqPercent="numeric",
                          seed="numeric",
                          mcmcIterations="numeric",
                          mcmcBurnin="numeric",
                          mcmcStep="numeric",
                          stopNames="ANY"),
           contains="GeneralSimulationsSummary")

################################################################################



##' Summarize the model-based design simulations, relative to a given truth
##'
##' @param object the \code{\linkS4class{Simulations}} object we want to
##' summarize
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity
##' @param target the target toxicity interval (default: 20-35\%) used for the
##' computations
##' @param \dots Additional arguments can be supplied here for \code{truth}
##' @return an object of class \code{\linkS4class{SimulationsSummary}}
##'
##' @example examples/Simulations-method-summary.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
            signature(object="Simulations"),
          def=
            function(object,
                     truth,
                     target=c(0.2, 0.35),
                     ...){
              
              ## call the parent method
              start <- callNextMethod(object=object,
                                      truth=truth,
                                      target=target,
                                      ...)
              
              doseGrid <- object@data[[1]]@doseGrid
              
              
              ## seed 
              seed<-object@seed
              
              ## mcmc options
              mcmcIterations <- object@mcmcIterations
              mcmcBurnin <- object@mcmcBurnin
              mcmcStep <- object@mcmcStep
              
              ## stopping rules class names
              stopNames<-object@stopNames
              
              
              
              ## dose level most often selected as MTD
              xMostSelected <-
                matchTolerance(start@doseMostSelected,
                               table=doseGrid)
              
              ## fitted toxicity rate at dose most often selected
              fitAtDoseMostSelected <-
                sapply(object@fit,
                       function(f){
                         f$middle[xMostSelected]
                       })
              
              ## mean fitted toxicity (average, lower and upper quantiles)
              ## at each dose level
              ## (this is required for plotting)
              meanFitMatrix <- sapply(object@fit,
                                      "[[",
                                      "middle")
              meanFit <- list(truth=
                                truth(doseGrid, ...),
                              average=
                                rowMeans(meanFitMatrix),
                              lower=
                                apply(meanFitMatrix,
                                      1L,
                                      quantile,
                                      0.025),
                              upper=
                                apply(meanFitMatrix,
                                      1L,
                                      quantile,
                                      0.975))
              
              ## median of model based median MTDs overall
              medianMTD <- median(object@MTD)
              
              ## median of model based median MTDs Rel. Error with regards to true MTD
              medianMTDRE <- median(object@MTDRE)
              
              ## mean of model CV of the MTDs (%) overall
              meanCVMTD <- mean(object@CV)
              
              ## TRUE MTD for scenario
              trueMTD <- object@trueMTD
              
              ## % trials (MTD >top dose) overall
              trialMTDaboveTopDose <- (sum((object@MTD >max(object@data[[1]]@doseGrid))=="TRUE")/length(object@data))*100
              
              ## % trials (MTD <first dose) overall
              trialMTDbelowFirstDose <- (sum((object@MTD <min(object@data[[1]]@doseGrid))=="TRUE")/length(object@data))*100
              
              ## number of cohorts per sim
              nCohorts <- sapply(object@data,function(d){max(d@cohort)})
              
              ## number of overdosed subjects with regards to the true MTD per sim
              nOD <- sapply(object@data,function(d){sum(d@x>object@trueMTD)})
              
              
              
              ## %trials each unique stopping rule fired individually 
              stopReasonsUnique <- sapply(object@stopReasonsSum,colSums)
              stopReasonsUniquePercent<- (rowSums(stopReasonsUnique)/length(stopReasonsUnique[1,]))*100
              
              ## %trials each stopping rule fired first based on the order of the stopping rule defined in the design  
              stopReasonsFirstHit<-apply(stopReasonsUnique,2,function(x){min(which(x == "1"))})
              stopReasonsFirstHitFreq<-table(stopReasonsFirstHit)
              
              if (dim(stopReasonsFirstHitFreq)<length(object@stopReasonsSum[[1]])){
                frame<-as.data.frame(stopReasonsFirstHitFreq)
                framefull<-data.frame(stopReasonsFirstHit=1:length(object@stopReasonsSum[[1]]))
                frame1 <- merge(x=framefull ,y= frame, by="stopReasonsFirstHit", all.x = TRUE)
                frame1$Freq[is.na(frame1$Freq)] = 0
                stopReasonsFirstHitFreq <- as.numeric(frame1$Freq)
              }
              ## percentage of different values of the freq 
              stopReasonsFirstHitFreqPercent <- (stopReasonsFirstHitFreq/sum(stopReasonsFirstHitFreq))*100
              
              
              ## give back an object of class SimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .SimulationsSummary(
                start,
                fitAtDoseMostSelected=fitAtDoseMostSelected,
                meanFit=meanFit,
                trialMTDaboveTopDose=trialMTDaboveTopDose,
                trialMTDbelowFirstDose=trialMTDbelowFirstDose,
                medianMTD=medianMTD,
                meanCVMTD=meanCVMTD,
                medianMTDRE=medianMTDRE,
                nCohorts=nCohorts,
                nOD=nOD,
                trueMTD=trueMTD,
                stopReasonsUniquePercent=stopReasonsUniquePercent,
                stopReasonsFirstHitFreqPercent=stopReasonsFirstHitFreqPercent,
                seed=seed,
                mcmcIterations=mcmcIterations,
                mcmcBurnin=mcmcBurnin,
                mcmcStep=mcmcStep,
                stopNames=stopNames
              )
              
              return(ret)
            })
