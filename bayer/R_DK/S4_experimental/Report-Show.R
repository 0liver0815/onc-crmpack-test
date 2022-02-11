##' A Reference Class to represent sequentially updated reporting objects.
##' @name Report
##' @field object The object from which to report
##' @field df the data frame to which columns are sequentially added
##' @field dfNames the names to which strings are sequentially added
Report <-
  setRefClass("Report",
              fields =
                list(object = "ANY",
                     df = "data.frame",
                     dfNames = "character"),
              methods = list(
                dfSave =
                  function(res, name) {
                    df <<- cbind(df, res)
                    dfNames <<- c(dfNames, name)
                    return(res)
                  },
                report =
                  function(slotName,
                           description,
                           percent=TRUE,
                           digits=1,
                           quantiles=c(0.1, 0.9),
                           subset=NULL,
                           doSum=FALSE,
                           method=2) {
                    vals <- slot(object, name=slotName)
                    if(! is.null(subset))
                      vals <- vals[subset,]  
                    if(doSum)
                      vals <- apply(vals, 2, sum)  
                    if(percent)
                    {
                      unit <- " %"
                      vals <- vals * 100
                    } else {
                      unit <- ""
                    }
                    
                    if (method==1){res <-
                      paste(round(mean(vals), digits), 
                            unit,
                            " (",
                            paste(round(quantile(vals,quantiles,na.rm=TRUE),
                                        digits),
                                  unit,
                                  collapse=", ",
                                  sep=""
                            ),
                            ")",
                            sep="")
                    }
                    if (method==2){res <-
                      paste(round(mean(vals), digits), 
                            unit,
                            " (",
                            paste(round(sd(vals,na.rm = TRUE),digits),
                                  unit
                            ),
                            "),",
                            " [",
                            paste(round(range(vals),digits),
                                  unit,
                                  collapse=", ",
                                  sep=""
                            ),
                            "]",
                            sep=""
                      )
                    }
                    
                    if (method==3){res <-
                      paste(round(mean(vals), digits), 
                            unit,
                            " (",
                            paste(round(sd(vals),digits),
                                  unit,
                            ),
                            "),",
                            " [",
                            paste(round(range(vals),digits),
                                  unit,
                                  collapse=", ",
                                  sep=""
                            ),
                            "]",
                            sep=""
                      )
                    }
                    
                    ## print result to the buffer
                    if (method==1){
                      cat(description, ":",
                          "mean",
                          dfSave(res, slotName),
                          "\n")
                    }
                    if (method==2){
                      cat(description, ":",
                          "mean (sd), [min, max]",
                          dfSave(res, slotName),
                          "\n")
                    }
                    if (method==3){
                      cat(description, ":",
                          dfSave(res, slotName),
                          "\n")
                    }
                  }))


##' Show the summary of the simulations
##'
##' @param object the \code{\linkS4class{GeneralSimulationsSummary}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @export
##' @keywords methods
setMethod("show",
          signature=
            signature(object="GeneralSimulationsSummary"),
          def=
            function(object){
              
              r <- Report$new(object=object,
                              df=
                                as.data.frame(matrix(nrow=1,
                                                     ncol=0)),
                              dfNames=character())
              
              cat("Summary of",
                  r$dfSave(object@nsim, "nsim"),
                  "simulations\n\n")
              
              cat("Random seed:",
                  r$dfSave(object@seed, "seed"),
                  "\n")
              
              cat("MCMC Iterations:",
                  r$dfSave(object@mcmcIterations, "mcmcIterations"),
                  "\n")
              cat("MCMC Burn-in:",
                  r$dfSave(object@mcmcBurnin, "mcmcBurnin"),
                  "\n")
              cat("MCMC Step:",
                  r$dfSave(object@mcmcStep, "mcmcStep"),
                  "\n\n")
              
              # cat("Target toxicity interval was",
              #     r$dfSave(paste(round(object@target * 100),
              #                    collapse=", "),
              #              "targetToxInterval"),
              #     "%\n")
              # cat("Target dose interval corresponding to this was",
              #     r$dfSave(paste(round(object@targetDoseInterval, 1),
              #                    collapse=", "),
              #              "targetDoseInterval"),
              #     "\n")
              # cat("Intervals are corresponding to",
              #     "10 and 90 % quantiles\n\n")
              # 
              # if(object@placebo){
              #   r$report("nObs",
              #            "Number of patients on placebo",
              #            percent=FALSE,
              #            subset=2)
              #   r$report("nObs",
              #            "Number of patients on active",
              #            percent=FALSE,
              #            subset=1)
              #   r$report("nObs",
              #            "Number of patients overall",
              #            percent=FALSE,
              #            doSum=TRUE)
              # }else{
              #   r$report("nObs",
              #            "Number of patients overall",
              #            percent=FALSE)
              # }
              # r$report("nAboveTarget",
              #          "Number of patients treated above target tox interval",
              #          percent=FALSE)
              # 
              # if(object@placebo){
              #   r$report("propDLTs",
              #            "Proportions of DLTs in the trials for patients on placebo",
              #            subset=2)
              #   r$report("propDLTs",
              #            "Proportions of DLTs in the trials for patients on active",
              #            subset=1)
              # }else{
              #   r$report("propDLTs",
              #            "Proportions of DLTs in the trials")  
              # }
              # r$report("meanToxRisk",
              #          "Mean toxicity risks for the patients on active")
              
              # r$report("doseSelected",
              #          "Doses selected as MTD",
              #          percent=FALSE, digits=1)
              
              # r$report("toxAtDosesSelected",
              #          "True toxicity at doses selected")
              
              # cat("Proportion of trials selecting target MTD:",
              #     r$dfSave(object@propAtTarget * 100,
              #              "percentAtTarget"),
              #     "%\n")
              
              # cat("Dose most often selected as MTD:",
              #     r$dfSave(object@doseMostSelected,
              #              "doseMostSelected"),
              #     "\n")
              
              # cat("Observed toxicity rate at dose most often selected:",
              #     r$dfSave(round(object@obsToxRateAtDoseMostSelected * 100),
              #              "obsToxRateAtDoseMostSelected"),
              #     "%\n")
              
              
              ## DK: get the individuals stop reasons
              stopLength<-length(object@stopReasonsUniquePercent)
              for (i in 1: stopLength){
                if(i<stopLength){
                  cat("Stops (individual) for reason",i,object@stopNames[i],":",
                      r$dfSave(round(object@stopReasonsUniquePercent[i], digits = 1),
                               "stopReasonsUniquePercent"),
                      "%\n")
                }else{
                  cat("Stops (individual) for reason",i, object@stopNames[i],":",
                      r$dfSave(round(object@stopReasonsUniquePercent[i], digits = 1),
                               "stopReasonsUniquePercent"),
                      "%\n\n")
                }
                
              }
              ## DK: get the first hit stop reasons
              stopLength<-length(object@stopReasonsFirstHitFreqPercent)
              for (i in 1: stopLength){
                if(i<stopLength){
                  cat("Stops (first hit) for reason",i,object@stopNames[i],":",
                      r$dfSave(round(object@stopReasonsFirstHitFreqPercent[i], digits = 1),
                               "stopReasonsFirstHitFreqPercent"),
                      "%\n")
                }else{
                  cat("Stops (first hit) for reason",i,object@stopNames[i],":",
                      r$dfSave(round(object@stopReasonsFirstHitFreqPercent[i], digits = 1),
                               "stopReasonsFirstHitFreqPercent"),
                      "%\n\n")
                }
                
              }
              
              cat("trials (MTD >top dose):",
                  r$dfSave(round(object@trialMTDaboveTopDose, digits = 1),
                           "trialMTDaboveTopDose"),
                  "%\n")
              
              cat("trials (MTD <first dose):",
                  r$dfSave(round(object@trialMTDbelowFirstDose, digits = 1),
                           "trialMTDbelowFirstDose"),
                  "%\n\n")
              
              cat("True MTD:",
                  r$dfSave(round(object@trueMTD, digits = 1),
                           "trueMTD"),
                  "\n")
              
              cat("Median MTD:",
                  r$dfSave(round(object@medianMTD, digits = 1),
                           "medianMTD"),
                  "\n")
              
              cat("Median MTD Rel. Error:",
                  r$dfSave(round(object@medianMTDRE, digits = 1),
                           "medianMTDRE"),
                  "%\n")
              
              cat("Mean CV (MTD):",
                  r$dfSave(round(object@meanCVMTD, digits = 1),
                           "meanCVMTD"),
                  "%\n\n")
              
              r$report("nOD",
                       "Number of over-dosed subjects",
                       percent=FALSE)
              
              if(object@placebo){
                r$report("nObs",
                         "Number of subjects on placebo",
                         percent=FALSE,
                         subset=2)
                r$report("nObs",
                         "Number of subjects on active",
                         percent=FALSE,
                         subset=1)
                r$report("nObs",
                         "Number of subjects overall",
                         percent=FALSE,
                         doSum=TRUE)
              }else{
                r$report("nObs",
                         "Number of subjects overall",
                         percent=FALSE)
              }
              
              r$report("nCohorts",
                       "Number of cohorts",
                       percent=FALSE)
              
              ## finally assign names to the df
              ## and return it invisibly
              names(r$df) <- r$dfNames
              invisible(r$df)
            })


##' Show the summary of the simulations
##'
##' @param object the \code{\linkS4class{SimulationsSummary}} object we want
##' to print
##' @return invisibly returns a data frame of the results with one row and
##' appropriate column names
##'
##' @example examples/Simulations-method-show-SimulationsSummary.R
##' @export
##' @keywords methods
setMethod("show",
          signature=
            signature(object="SimulationsSummary"),
          def=
            function(object){
              
              ## call the parent method
              df <- callNextMethod(object)
              dfNames <- names(df)
              
              ## start report object
              r <- Report$new(object=object,
                              df=df,
                              dfNames=dfNames)
              
              # ## add one reporting line
              # r$report("fitAtDoseMostSelected",
              #          "Fitted toxicity rate at dose most often selected")
              # 
              # 
              
              ## and return the updated information
              names(r$df) <- r$dfNames
              invisible(r$df)
            })
