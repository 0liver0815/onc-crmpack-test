
#' @describeIn examine Examine a model-based Bayer CRM
#'
#' @param mcmcOptions object of class \code{\linkS4class{McmcOptions}},
#' giving the MCMC options for each evaluation in the trial. By default,
#' the standard options are used
#' 
#' @example examples/design-method-examine-Design.R
setMethod("examine",
          signature=
            signature(object="BayDesign"),
          def=
            function(object, 
                     mcmcOptions=McmcOptions(), 
                     ...,
                     maxNoIncrement){
              
              # start with the empty table
              ret <- data.frame(dose=numeric(),
                                DLTs=integer(),
                                nextDose=numeric(),
                                stop=logical(),
                                increment=integer())
              
              # start the base data with the provided one
              baseData <- object@data
              
              # are we finished and can stop?
              stopit <- FALSE
              
              # counter how many contiguous doses at 0 DLTs with 
              # no increment
              noIncrementCounter <- 0L
              
              # what is the next dose to be used?
              # initialize with starting dose
              thisDose <- object@startingDose
              
              # inside this loop we continue filling up the table, until
              # stopping
              while(! stopit)
              {
                # what is the cohort size at this dose?
                thisSize <- size(cohortSize=object@cohortSize,
                                 dose=thisDose,
                                 data=baseData)
                
                if(baseData@placebo)
                  thisSize.PL <- size(cohortSize=object@PLcohortSize,
                                      dose=thisDose,
                                      data=baseData)
                
                # for all possible number of DLTs:
                for(numDLTs in 0:thisSize)
                {
                  # update data with corresponding DLT vector
                  if(baseData@placebo && (thisSize.PL > 0L) ){
                    thisData <- update(object=baseData,
                                       x=baseData@doseGrid[1],
                                       y=rep(0,thisSize.PL))
                    
                    thisData <-
                      update(object=thisData,
                             x=thisDose,
                             y=
                               rep(x=c(0, 1),
                                   times=
                                     c(thisSize - numDLTs,
                                       numDLTs)),
                             newCohort=FALSE)
                    
                  }else{
                    thisData <-
                      update(object=baseData,
                             x=thisDose,
                             y=
                               rep(x=c(0, 1),
                                   times=
                                     c(thisSize - numDLTs,
                                       numDLTs)))
                  }
                  
                  # what is the dose limit?
                  doselimit <- maxDose(object@increments,
                                       data=thisData)
                  
                  # generate samples from the model
                  thisSamples <- mcmc(data=thisData,
                                      model=object@model,
                                      options=mcmcOptions)
                  
                  # => what is the next best dose?
                  nextDose <- nextBest(object@nextBest,
                                       doselimit=doselimit,
                                       samples=thisSamples,
                                       model=object@model,
                                       data=thisData)$value
                  
                  # compute relative increment in percent
                  thisIncrement <-
                    round((nextDose - thisDose) / thisDose * 100)
                  
                  # evaluate stopping rules
                  stopThisTrial <- stopTrial(object@stopping,
                                             dose=nextDose,
                                             samples=thisSamples,
                                             model=object@model,
                                             data=thisData)
                  
                  # append information to the data frame
                  ret <- rbind(ret,
                               list(dose=thisDose,
                                    DLTs=numDLTs,
                                    nextDose=nextDose,
                                    stop=stopThisTrial,
                                    increment=as.integer(thisIncrement)))
                }
                
                # change base data
                if(baseData@placebo && (thisSize.PL > 0L) ){
                  baseData <-
                    update(object=baseData,
                           x=baseData@doseGrid[1],
                           y=rep(0, thisSize.PL))
                  
                  baseData <-
                    update(object=baseData,
                           x=thisDose,
                           y=rep(0, thisSize),
                           newCohort=FALSE)
                  
                }else{
                  baseData <-
                    update(object=baseData,
                           x=thisDose,
                           y=rep(0, thisSize))
                }
                
                # what are the results if 0 DLTs?
                resultsNoDLTs <- subset(tail(ret, thisSize + 1),
                                        dose==thisDose & DLTs==0)
                
                # what is the new dose then accordingly?
                newDose <- as.numeric(resultsNoDLTs$nextDose)
                
                # what is the difference to the previous dose?
                doseDiff <- newDose - thisDose
                
                # update the counter for no increments of the dose
                if(doseDiff == 0)
                {
                  noIncrementCounter <- noIncrementCounter + 1L
                } else {
                  noIncrementCounter <- 0L
                }
                
                # would stopping rule be fulfilled already?
                stopAlready <- resultsNoDLTs$stop
                
                # update dose
                thisDose <- newDose
                
                # too many times no increment?
                stopNoIncrement <- (noIncrementCounter >= maxNoIncrement)
                if(stopNoIncrement) 
                  warning(paste("Stopping because", 
                                noIncrementCounter,
                                "times no increment vs. previous dose"))
                
                # check if we can stop:
                # either when we have reached the highest dose in the
                # next cohort, or when the stopping rule is already 
                # fulfilled, or when too many times no increment
                stopit <- (thisDose >= max(object@data@doseGrid)) ||
                  stopAlready || stopNoIncrement
              }
              
              return(ret)
            })
