#' Simulate outcomes from a Bayer CRM design
#'
#' @param object the \code{\linkS4class{BayDesign}} object we want to simulate
#' data from
#' @param nsim the number of simulations (default: 1)
#' @param seed see \code{\link{setSeed}}
#' @param truth a function which takes as input a dose (vector) and returns the
#' true probability (vector) for toxicity. Additional arguments can be supplied
#' in \code{args}.
#' @param args data frame with arguments for the \code{truth} function. The
#' column names correspond to the argument names, the rows to the values of the
#' arguments. The rows are appropriately recycled in the \code{nsim}
#' simulations. In order to produce outcomes from the posterior predictive
#' distribution, e.g, pass an \code{object} that contains the data observed so
#' far, \code{truth} contains the \code{prob} function from the model in
#' \code{object}, and \code{args} contains posterior samples from the model.
#' @param firstSeparate enroll the first patient separately from the rest of
#' the cohort? (not default) If yes, the cohort will be closed if a DLT occurs
#' in this patient.
#' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
#' giving the MCMC options for each evaluation in the trial. By default,
#' the standard options are used
#' @param parallel should the simulation runs be parallelized across the
#' clusters of the computer? (not default)
#' @param nCores how many cores should be used for parallel computing?
#' Defaults to the number of cores on the machine, maximum 5.
#' @param \dots not used
#'
#' @return an object of class \code{\linkS4class{Simulations}}
#'
#' @example examples/design-method-simulate-Design.R
#' @export
#' @importFrom parallel detectCores
#' @keywords methods
setMethod("simulate",
          signature=
            signature(object="BayDesign",
                      nsim="ANY",
                      seed="ANY"),
          def=
            function(object, nsim=1L, seed=NULL,
                     truth, args=NULL, firstSeparate=FALSE,
                     mcmcOptions=McmcOptions(),
                     parallel=FALSE, nCores=
                       min(parallel::detectCores(), 32),
                     ...){
              
              nsim <- safeInteger(nsim)
              
              # checks and extracts
              stopifnot(is.function(truth),
                        is.bool(firstSeparate),
                        is.scalar(nsim),
                        nsim > 0,
                        is.bool(parallel),
                        is.scalar(nCores),
                        nCores > 0)
              
              args <- as.data.frame(args)
              nArgs <- max(nrow(args), 1L)
              
              # seed handling
              RNGstate <- setSeed(seed)
              
              # from this,
              # generate the individual seeds for the simulation runs
              simSeeds <- sample(x=seq_len(1e5), size=nsim)
              
              # the function to produce the run a single simulation
              # with index "iterSim"
              runSim <- function(iterSim)
              {
                # set the seed for this run
                set.seed(simSeeds[iterSim])
                
                # what is now the argument for the truth?
                # (appropriately recycled)
                thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop=FALSE]
                
                # so this truth is...
                thisTruth <- function(dose)
                {
                  do.call(truth,
                          # First argument: the dose
                          c(dose,
                            # Following arguments
                            thisArgs))
                }
                
                # start the simulated data with the provided one
                thisData <- object@data
                
                # In case there are placebo
                if(thisData@placebo)
                  # what is the probability for tox. at placebo?
                  thisProb.PL <- thisTruth(object@data@doseGrid[1])
                
                # shall we stop the trial?
                # First, we want to continue with the starting dose.
                # This variable is updated after each cohort in the loop.
                stopit <- FALSE
                
                # what is the next dose to be used?
                # initialize with starting dose
                thisDose <- object@startingDose
                
                # inside this loop we simulate the whole trial, until stopping
                while(! stopit)
                {
                  # what is the probability for tox. at this dose?
                  thisProb <- thisTruth(thisDose)
                  
                  # what is the cohort size at this dose?
                  thisSize <- size(cohortSize=object@cohortSize,
                                   dose=thisDose,
                                   data=thisData)
                  
                  # In case there are placebo
                  if(thisData@placebo)
                    thisSize.PL <- size(cohortSize=object@PLcohortSize,
                                        dose=thisDose,
                                        data=thisData)
                  
                  
                  # simulate DLTs: depends on whether we
                  # separate the first patient or not.
                  if(firstSeparate && (thisSize > 1L))
                  {
                    # dose the first patient
                    thisDLTs <- rbinom(n=1L,
                                       size=1L,
                                       prob=thisProb)
                    
                    if( thisData@placebo && (thisSize.PL > 0L) )
                      thisDLTs.PL <- rbinom(n=1L,
                                            size=1L,
                                            prob=thisProb.PL)
                    
                    # if there is no DLT:
                    if(thisDLTs == 0)
                    {
                      # enroll the remaining patients
                      thisDLTs <- c(thisDLTs,
                                    rbinom(n=thisSize - 1L,
                                           size=1L,
                                           prob=thisProb))
                      
                      if( thisData@placebo && (thisSize.PL > 0L) )
                        thisDLTs.PL <- c(thisDLTs.PL,
                                         rbinom(n=thisSize.PL,
                                                size=1L,
                                                prob=thisProb.PL))
                    }
                    
                  } else {
                    # we can directly dose all patients
                    thisDLTs <- rbinom(n=thisSize,
                                       size=1L,
                                       prob=thisProb)
                    
                    if( thisData@placebo && (thisSize.PL > 0L) )
                      thisDLTs.PL <- rbinom(n=thisSize.PL,
                                            size=1L,
                                            prob=thisProb.PL) 
                  }
                  
                  # update the data with this placebo (if any) cohort and then with active dose
                  if( thisData@placebo && (thisSize.PL > 0L) ){
                    thisData <- update(object=thisData,
                                       x=object@data@doseGrid[1],
                                       y=thisDLTs.PL)
                    
                    # update the data with active dose
                    thisData <- update(object=thisData,
                                       x=thisDose,
                                       y=thisDLTs,
                                       newCohort=FALSE)
                  }else{
                    # update the data with this cohort
                    thisData <- update(object=thisData,
                                       x=thisDose,
                                       y=thisDLTs)
                  }
                  
                  # what is the dose limit?
                  doselimit <- maxDose(object@increments,
                                       data=thisData)
                  
                  # generate samples from the model
                  thisSamples <- mcmc(data=thisData,
                                      model=object@model,
                                      options=mcmcOptions)
                  
                  # => what is the next best dose?
                  thisDose <- nextBest(object@nextBest,
                                       doselimit=doselimit,
                                       samples=thisSamples,
                                       model=object@model,
                                       data=thisData)$value
                  
                  # evaluate stopping rules
                  stopit <- stopTrial(object@stopping,
                                      dose=thisDose,
                                      samples=thisSamples,
                                      model=object@model,
                                      data=thisData)
                }
                
                # get the fit
                thisFit <- fit(object=thisSamples,
                               model=object@model,
                               data=thisData)
                
                # get the allocation criteria at the end of the trial
                thisAllocation <- nextBest(object@nextBest,
                                           doselimit=doselimit,
                                           samples=thisSamples,
                                           model=object@model,
                                           data=thisData)$allocation
                
                # get the MTD estimate from the samples
                baythisMTDsamp <- dose(object@nextBest@target,
                                       model=object@model,
                                       samples = thisSamples)
                
                #create a list with the estimates for the MTD and allocation citeria
                baythisMTD <- list(medianMTD=median(baythisMTDsamp),
                                   CVMTD=mad(baythisMTDsamp)/median(baythisMTDsamp),
                                   allocation=thisAllocation)
                
                #add the MTD to the fit in the last column
                #thisFit <- (rbind(thisFit, baythisMTD))
                
                # return the results
                thisResult <-
                  list(data=thisData,
                       dose=thisDose,
                       fit=
                         subset(thisFit,
                                select=c(middle, lower, upper)),
                       estimates=baythisMTD,
                       stop=
                         attr(stopit,
                              "message"))
                return(thisResult)
              }
              
              resultList <- getResultList(fun=runSim,
                                          nsim=nsim,
                                          vars=
                                            c("simSeeds",
                                              "args",
                                              "nArgs",
                                              "firstSeparate",
                                              "truth",
                                              "object",
                                              "mcmcOptions"),
                                          parallel=if(parallel) nCores else NULL)
              
              # put everything in the Simulations format:
              
              # setup the list for the simulated data objects
              dataList <- lapply(resultList, "[[", "data")
              
              # the vector of the final dose recommendations
              recommendedDoses <- as.numeric(sapply(resultList, "[[", "dose"))
              
              # setup the list for the final fits
              fitList <- lapply(resultList, "[[", "fit")
              
              # setup the list for the final estimates
              estimatesList <- lapply(resultList, "[[", "estimates")
              
              # the reasons for stopping
              stopReasons <- lapply(resultList, "[[", "stop")
              
              # return the results in the Simulations class object
              ret <- SimulationsBay(data=dataList,
                                    doses=recommendedDoses,
                                    fit=fitList,
                                    estimates=estimatesList,
                                    stopReasons=stopReasons,
                                    seed=RNGstate)
              
              return(ret)
            })


