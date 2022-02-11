##' General class for the simulations output Bayer
##'
##' This class captures trial simulations.
##'
##' Here also the random generator state before starting the simulation is
##' saved, in order to be able to reproduce the outcome. For this just use
##' \code{\link{set.seed}} with the \code{seed} as argument before running
##' \code{\link{simulate,Design-method}}.
##'
##' @slot data list of produced \code{\linkS4class{Data}} objects
##' @slot doses the vector of final dose recommendations
##' @slot seed random generator state before starting the simulation
##'
##' @export
##' @keywords classes
.GeneralSimulationsBay <-
  setClass(Class="GeneralSimulationsBay",
           representation(data="list",
                          doses="numeric",
                          seed="integer"),
           prototype(data=
                       list(Data(x=1:2,
                                 y=0:1,
                                 doseGrid=1:2,
                                 ID=1:2,
                                 cohort=1:2),
                            Data(x=3:4,
                                 y=0:1,
                                 doseGrid=3:4,
                                 ID=3:4,
                                 cohort=3:4)),
                     doses=c(1, 2),
                     seed=1L),
           validity=
             function(object){
               o <- Validate()
               
               nSims <- length(object@data)
               
               o$check(all(sapply(object@data, is, "Data")),
                       "all data elements must be Data objects")
               o$check(identical(length(object@doses), nSims),
                       "doses must have same length as the data list")
               
               o$result()
             })
validObject(.GeneralSimulationsBay())


##' Initialization function for "GeneralSimulations"
##'
##' @param data see \code{\linkS4class{GeneralSimulations}}
##' @param doses see \code{\linkS4class{GeneralSimulations}}
##' @param seed see \code{\linkS4class{GeneralSimulations}}
##' @return the \code{\linkS4class{GeneralSimulations}} object
##'
##' @export
##' @keywords methods
GeneralSimulationsBay <- function(data,
                                  doses,
                                  seed)
{
  .GeneralSimulationsBay(data=data,
                         doses=doses,
                         seed=safeInteger(seed))
}




##' Class for the summary of general simulations output
##'
##' Note that objects should not be created by users, therefore no
##' initialization function is provided for this class.
##'
##' @slot target target toxicity interval
##' @slot targetDose corresponding target dose interval
##' @slot targetDoseLevel corresponding target dose interval
##' @slot nsim number of simulations
##' @slot propDLTs proportions of DLTs in the trials
##' @slot meanToxRisk mean toxicity risks for the patients
##' @slot doseSelected doses selected as MTD
##' @slot toxAtDosesSelected true toxicity at doses selected
##' @slot propAtTarget Proportion of trials selecting target MTD
##' @slot doseMostSelected dose most often selected as MTD
##' @slot obsToxRateAtDoseMostSelected observed toxicity rate at dose most often
##' selected
##' @slot nObs number of patients overall
##' @slot nAboveTarget number of patients treated above target tox interval
##' @slot doseGrid the dose grid that has been used
##' @slot placebo set to TRUE (default is FALSE) for a design with placebo
##'
##' @export
##' @keywords classes
.GeneralSimulationsSummaryBay <-
  setClass(Class="GeneralSimulationsSummaryBay",
           representation(target="numeric",
                          targetDose="character",
                          targetDoseLevel="numeric",
                          nsim="integer",
                          propDLTs="ANY",
                          meanToxRisk="numeric",
                          doseSelected="numeric",
                          toxAtDosesSelected="numeric",
                          propAtTarget="numeric",
                          doseMostSelected="numeric",
                          obsToxRateAtDoseMostSelected="numeric",
                          nObs="ANY",
                          nAboveTarget="integer",
                          doseGrid="numeric",
                          placebo="logical"))


##' Class for the summary of model-based simulations output
##'
##' In addition to the slots in the parent class
##' \code{\linkS4class{GeneralSimulationsSummaryBay}}, it contains two slots with
##' model fit information.
##'
##' Note that objects should not be created by users, therefore no
##' initialization function is provided for this class.
##'
##' @slot fitAtDoseMostSelected fitted toxicity rate at dose most often selected
##' @slot meanFit list with the average, lower (2.5%) and upper (97.5%)
##' quantiles of the mean fitted toxicity at each dose level
##' @slot medianMTD list with the medianMTD
##' @slot cvMTD list with the CV of the MTD
##' @slot stopreason stop reason triggered
##' @slot stopreason.unique only one stop reason counted
##'
##' @export
##' @keywords classes
.SimulationsSummaryBay <-
  setClass(Class="SimulationsSummaryBay",
           representation(fitAtDoseMostSelected="numeric",
                          meanFit="list",
                          medianMTD="numeric",
                          cvMTD="numeric",
                          stopreason="numeric",
                          stopreason.unique="table"),
           contains="GeneralSimulationsSummaryBay")





##' Summarize the simulations, relative to a given truth
##'
##' @param object the \code{\linkS4class{GeneralSimulations}} object we want to
##' summarize
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity
##' @param target the target toxicity interval (default: 20-35\%) used for the
##' computations
##' @param \dots Additional arguments can be supplied here for \code{truth}
##' @return an object of class \code{\linkS4class{GeneralSimulationsSummary}}
##'
##' @export
##' @keywords methods
setMethod("summary",
          signature=
            signature(object="GeneralSimulationsBay"),
          def=
            function(object,
                     truth,
                     target=0.3,
                     ...){
              
              ## extract dose grid
              doseGrid <- object@data[[1]]@doseGrid
              
              ## evaluate true toxicity at doseGrid
              trueTox <- truth(doseGrid, ...)
              
              ## find dose level and dose that is closest to target
              targetDoseLevel <- max(which(abs(trueTox-target)==min(abs(trueTox-target))))
              if (trueTox[1] > target){
                targetDoseLevel <- 0
                targetDose <- paste("<", doseGrid[1])
              }else if (tail(trueTox,1) < target){
                targetDoseLevel <- length(doseGrid)+1
                targetDose <- paste(">", tail(doseGrid,1))
              }else {targetDose <- paste(doseGrid[targetDoseLevel])} 
              
              ## what are the levels above target interval?
              xAboveTarget <- which(trueTox > target)
              
              ## proportion of DLTs in a trial:
              if(object@data[[1]]@placebo){
                if( sum(object@data[[1]]@x == doseGrid[1]) ){
                  propDLTs <- sapply(object@data,
                                     function(d){
                                       tapply(d@y,
                                              factor(d@x == d@doseGrid[1], 
                                                     labels=c("ACTV","PLCB")), 
                                              mean)
                                     })
                }else{
                  propDLTs <- sapply(object@data,
                                     function(d){
                                       c('ACTV' = mean(d@y),'PLCB' = NA)   
                                     }) 
                }
              }else{
                propDLTs <- sapply(object@data,
                                   function(d){
                                     mean(d@y)   
                                   })
              }
              
              ## mean toxicity risk
              if(object@data[[1]]@placebo){
                meanToxRisk <- sapply(object@data,
                                      function(d){
                                        mean(trueTox[d@xLevel[d@xLevel != 1]])
                                      })
              }else{
                meanToxRisk <- sapply(object@data,
                                      function(d){
                                        mean(trueTox[d@xLevel])
                                      })
              }
              
              ## doses selected for MTD
              doseSelected <- object@doses
              
              ## replace NA by 0
              doseSelected[is.na(doseSelected)] <- 0
              
              ## dose most often selected as MTD
              doseMostSelected <-
                as.numeric(names(which.max(table(doseSelected))))
              xMostSelected <-
                matchTolerance(doseMostSelected,
                               table=doseGrid)
              
              ## observed toxicity rate at dose most often selected
              ## Note: this does not seem very useful!
              ## Reason: In case of a fine grid, few patients if any
              ## will have been treated at this dose.
              tmp <-
                sapply(object@data,
                       function(d){
                         whichAtThisDose <- which(d@x == doseMostSelected)
                         nAtThisDose <- length(whichAtThisDose)
                         nDLTatThisDose <- sum(d@y[whichAtThisDose])
                         return(c(nAtThisDose=nAtThisDose,
                                  nDLTatThisDose=nDLTatThisDose))
                       })
              
              obsToxRateAtDoseMostSelected <-
                mean(tmp["nDLTatThisDose",]) / mean(tmp["nAtThisDose",])
              
              ## number of patients overall
              if(object@data[[1]]@placebo){
                nObs <- sapply(object@data,
                               function(x){
                                 data.frame(n.ACTV = sum(x@xLevel != 1L),
                                            n.PLCB = sum(x@xLevel == 1L))
                               })
                nObs <- matrix(unlist(nObs), dim(nObs))
              }else{
                nObs <- sapply(object@data,
                               slot,
                               "nObs")
              }
              
              
              ## number of patients treated above target tox interval
              nAboveTarget <- sapply(object@data,
                                     function(d){
                                       sum(d@xLevel %in% xAboveTarget)
                                     })
              
              ## Tox at dose level selected
              toxAtDoses <- truth(doseSelected, ...)
              #propAtTarget <- mean((toxAtDoses > target[1]) &
              #                       (toxAtDoses < target[2]))
              ## Proportion of trials selecting target MTD
              propAtTarget <- mean(doseSelected == doseGrid[targetDoseLevel])
              
              ## give back an object of class GeneralSimulationsSummary,
              ## for which we then define a print / plot method
              ret <-
                .GeneralSimulationsSummaryBay(
                  target=target,
                  targetDoseLevel=targetDoseLevel,
                  targetDose=targetDose,
                  nsim=length(object@data),
                  propDLTs=propDLTs,
                  meanToxRisk=meanToxRisk,
                  doseSelected=doseSelected,
                  doseMostSelected=doseMostSelected,
                  obsToxRateAtDoseMostSelected=obsToxRateAtDoseMostSelected,
                  nObs=nObs,
                  nAboveTarget=nAboveTarget,
                  toxAtDosesSelected=toxAtDoses,
                  propAtTarget=propAtTarget,
                  doseGrid=doseGrid,
                  placebo=object@data[[1]]@placebo)  
              
              
              return(ret)
            })

