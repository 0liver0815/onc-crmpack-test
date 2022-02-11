##' Show the summary of the simulations
##'
##' @param object the \code{\linkS4class{GeneralSimulationsSummaryBay}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @export
##' @keywords methods
setMethod("show",
          signature=
            signature(object="GeneralSimulationsSummaryBay"),
          def=
            function(object){
              
              r <- ReportBay$new(object=object,
                                 df=
                                   as.data.frame(matrix(nrow=1,
                                                        ncol=0)),
                                 dfNames=character())
              
              cat("Summary of",
                  r$dfSave(object@nsim, "nsim"),
                  "simulations\n\n")
              
              cat("Target toxicity was",
                  r$dfSave(paste(round(object@target * 100),
                                 collapse=" "),
                           "targetTox"),
                  "%\n")
              cat("Target dose level corresponding to this was",
                  r$dfSave(paste(round(object@targetDoseLevel, 1),
                                 collapse=" "),
                           "targetDoseLevel"),
                  "\n")
              cat("Target dose corresponding to this was",
                  r$dfSave(object@targetDose,
                           "targetDose"),
                  "\n")
              cat("Intervals are corresponding to",
                  "0, 10, 90 and 100 % quantiles\n\n")
              
              if(object@placebo){
                r$report("nObs",
                         "Number of patients on placebo",
                         percent=FALSE,
                         subset=2)
                r$report("nObs",
                         "Number of patients on active",
                         percent=FALSE,
                         subset=1)
                r$report("nObs",
                         "Number of patients overall",
                         percent=FALSE,
                         doSum=TRUE)
              }else{
                r$report("nObs",
                         "Number of patients overall",
                         percent=FALSE)
              }
              r$report("nAboveTarget",
                       "Number of patients treated above target tox",
                       percent=FALSE)
              
              if(object@placebo){
                r$report("propDLTs",
                         "Proportions of DLTs in the trials for patients on placebo",
                         subset=2)
                r$report("propDLTs",
                         "Proportions of DLTs in the trials for patients on active",
                         subset=1)
              }else{
                r$report("propDLTs",
                         "Proportions of DLTs in the trials")  
              }
              r$report("meanToxRisk",
                       "Mean toxicity risks for the patients on active")
              r$report("doseSelected",
                       "Doses selected as MTD",
                       percent=FALSE, digits=1)
              r$report("toxAtDosesSelected",
                       "True toxicity at doses selected")
              cat("Proportion of trials selecting target MTD:",
                  r$dfSave(object@propAtTarget * 100,
                           "percentAtTarget"),
                  "%\n")
              cat("Dose most often selected as MTD:",
                  r$dfSave(object@doseMostSelected,
                           "doseMostSelected"),
                  "\n")
              cat("Observed toxicity rate at dose most often selected:",
                  r$dfSave(round(object@obsToxRateAtDoseMostSelected * 100),
                           "obsToxRateAtDoseMostSelected"),
                  "%\n")
              
              ## finally assign names to the df
              ## and return it invisibly
              names(r$df) <- r$dfNames
              invisible(r$df)
            })
