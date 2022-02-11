
################################################################################
getResultList <- function(fun,
                          nsim,
                          vars,
                          parallel=NULL)
{
  ret <-
    if(is.null(parallel))
    {
      lapply(X=seq_len(nsim),
             FUN=fun)
    } else {
      
      ## check that parallel parameter makes sense
      stopifnot(is.scalar(parallel), parallel > 0)
      
      ## now process all simulations
      cores <- min(safeInteger(parallel),
                   min(parallel::detectCores(), 5))
      
      ## start the cluster
      cl <- parallel::makeCluster(cores)
      
      ## load the required R package
      parallel::clusterEvalQ(cl, {
        library(crmPack)
        NULL
      })
      
      ## export local variables
      parallel::clusterExport(cl=cl,
                              varlist=vars,
                              envir=parent.frame())
      ## parent.frame() gives back the caller environment
      ## (different from parent.env() which returns
      ## the environment where this function has been
      ## defined!)
      
      ## export all global variables
      parallel::clusterExport(cl=cl,
                              varlist=ls(.GlobalEnv))
      
      ## now do the computations
      res <- parallel::parLapply(cl=cl,
                                 X=seq_len(nsim),
                                 fun=fun)
      
      ## stop the cluster
      parallel::stopCluster(cl)
      
      res
    }
  
  return(ret)
}

################################################################################


################################################################################
.Simulations <-
  setClass(Class="Simulations",
           representation(fit="list",
                          stopReasons="list",
                          stopReasonsSum="list",
                          MTD="numeric",
                          CV="numeric",
                          MTDRE="numeric",
                          trueMTD="numeric",
                          mcmcIterations="numeric",
                          mcmcBurnin="numeric",
                          mcmcStep="numeric",
                          stopNames="ANY",
                          seed_iterSim="numeric"
           ),
           ## note: this prototype is put together with the prototype
           ## for GeneralSimulations
           prototype(fit=
                       list(c(0.1, 0.2),
                            c(0.1, 0.2)),
                     stopReasons=
                       list("A", "A"),
                     stopReasonsSum=
                       list("FALSE", "TRUE"),
                     MTD=c(1, 2),
                     CV=c(1, 2),
                     MTDRE=c(1, 2),
                     trueMTD=2,
                     mcmcIterations=2000,
                     mcmcBurnin=1000,
                     mcmcStep=1,
                     stopNames=list("A", "A"),
                     seed_iterSim=c(1, 2)
           ),
           contains="GeneralSimulations",
           validity=
             function(object){
               o <- Validate()
               
               nSims <- length(object@data)
               
               o$check(identical(length(object@fit), nSims),
                       "fit must have same length as data")
               o$check(identical(length(object@stopReasons), nSims),
                       "stopReasons must have same length as data")
               o$check(identical(length(object@stopReasonsSum), nSims),
                       "stopReasonsSum must have same length as data")
               o$check(identical(length(object@MTD), nSims),
                       "MTDs must have same length as the data")
               o$check(identical(length(object@CV), nSims),
                       "CVs must have same length as the data")
               o$check(identical(length(object@MTDRE), nSims),
                       "median MTD relative errors must have same length as the data")
               
               # o$check(identical(length(object@seed_iterSim), nSims),
               #         "seed_iterSim must have same length as the data")
               
               
               o$result()
             })
validObject(.Simulations())



Simulations <- function(fit,
                        stopReasons,stopReasonsSum,MTD,CV,MTDRE,trueMTD,mcmcIterations,mcmcBurnin,mcmcStep,stopNames,seed_iterSim,
                        ...)
{
  start <- GeneralSimulations(...)
  .Simulations(start,
               fit=fit,
               stopReasons=stopReasons,
               stopReasonsSum=stopReasonsSum,
               MTD=MTD,
               CV=CV,
               MTDRE=MTDRE,
               trueMTD=trueMTD,
               mcmcIterations=mcmcIterations,
               mcmcBurnin=mcmcBurnin,
               mcmcStep=mcmcStep,
               stopNames=stopNames,
               seed_iterSim=seed_iterSim
  )
}



################################################################################

setMethod("simulate",
          signature=
            signature(object="Design",
                      nsim="ANY",
                      seed="ANY"),
          def=
            function(object, nsim=1L, seed=NULL,
                     truth, trueMTD=NULL, args=NULL, firstSeparate=FALSE,
                     mcmcOptions=McmcOptions(),
                     parallel=FALSE, nCores=
                       min(parallel::detectCores(), 5),
                     ...){
              
              nsim <- safeInteger(nsim)
              
              ## checks and extracts
              stopifnot(is.function(truth),
                        is.bool(firstSeparate),
                        is.scalar(nsim),
                        nsim > 0,
                        is.bool(parallel),
                        is.scalar(nCores),
                        nCores > 0)
              
              args <- as.data.frame(args)
              nArgs <- max(nrow(args), 1L)
              
              ## seed handling
              RNGstate <- setSeed(seed)
              
              ## mcmc options
              mcmcIterations<- mcmcOptions@iterations
              mcmcBurnin<-mcmcOptions@burnin
              mcmcStep<- mcmcOptions@step
              
              ## stopping rules class names
              stopLength<-length(object@stopping@stopList)
              stopNames <- rep(NA, stopLength)
              
              for (i in 1: stopLength){
                stopNames[i]<- class(object@stopping@stopList[[i]])[1]
                
              }
              
              ## from this,
              ## generate the individual seeds for the simulation runs
              simSeeds <- sample(x=seq_len(1e5), size=nsim)
              
              ## the function to produce the run a single simulation
              ## with index "iterSim"
              runSim <- function(iterSim)
              {
                ## set the seed for this run
                set.seed(simSeeds[iterSim])
                
                
                ## what is now the argument for the truth?
                ## (appropriately recycled)
                thisArgs <- args[(iterSim - 1) %% nArgs + 1, , drop=FALSE]
                
                ## so this truth is...
                thisTruth <- function(dose)
                {
                  do.call(truth,
                          ## First argument: the dose
                          c(dose,
                            ## Following arguments
                            thisArgs))
                }
                
                ## start the simulated data with the provided one
                thisData <- object@data
                
                # In case there are placebo
                if(thisData@placebo)
                  ## what is the probability for tox. at placebo?
                  thisProb.PL <- thisTruth(object@data@doseGrid[1])
                
                ## shall we stop the trial?
                ## First, we want to continue with the starting dose.
                ## This variable is updated after each cohort in the loop.
                stopit <- FALSE
                
                ## what is the next dose to be used?
                ## initialize with starting dose
                thisDose <- object@startingDose
                
                ## inside this loop we simulate the whole trial, until stopping
                while(! stopit)
                {
                  ## what is the probability for tox. at this dose?
                  thisProb <- thisTruth(thisDose)
                  
                  ## what is the cohort size at this dose?
                  thisSize <- size(cohortSize=object@cohortSize,
                                   dose=thisDose,
                                   data=thisData)
                  
                  ## In case there are placebo
                  if(thisData@placebo)
                    thisSize.PL <- size(cohortSize=object@PLcohortSize,
                                        dose=thisDose,
                                        data=thisData)
                  
                  
                  ## simulate DLTs: depends on whether we
                  ## separate the first patient or not.
                  if(firstSeparate && (thisSize > 1L))
                  {
                    ## dose the first patient
                    thisDLTs <- rbinom(n=1L,
                                       size=1L,
                                       prob=thisProb)
                    
                    if( thisData@placebo && (thisSize.PL > 0L) )
                      thisDLTs.PL <- rbinom(n=1L,
                                            size=1L,
                                            prob=thisProb.PL)
                    
                    ## if there is no DLT:
                    if(thisDLTs == 0)
                    {
                      ## enroll the remaining patients
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
                    ## we can directly dose all patients
                    thisDLTs <- rbinom(n=thisSize,
                                       size=1L,
                                       prob=thisProb)
                    
                    if( thisData@placebo && (thisSize.PL > 0L) )
                      thisDLTs.PL <- rbinom(n=thisSize.PL,
                                            size=1L,
                                            prob=thisProb.PL) 
                  }
                  
                  ## update the data with this placebo (if any) cohort and then with active dose
                  if( thisData@placebo && (thisSize.PL > 0L) ){
                    thisData <- update(object=thisData,
                                       x=object@data@doseGrid[1],
                                       y=thisDLTs.PL)
                    
                    ## update the data with active dose
                    thisData <- update(object=thisData,
                                       x=thisDose,
                                       y=thisDLTs,
                                       newCohort=FALSE)
                  }else{
                    ## update the data with this cohort
                    thisData <- update(object=thisData,
                                       x=thisDose,
                                       y=thisDLTs)
                  }
                  
                  ## what is the dose limit?
                  doselimit <- maxDose(object@increments,
                                       data=thisData)
                  
                  ## generate samples from the model
                  thisSamples <- mcmc(data=thisData,
                                      model=object@model,
                                      options=mcmcOptions)
                  
                  ## => what is the next best dose?
                  thisDose <- nextBest(object@nextBest,
                                       doselimit=doselimit,
                                       samples=thisSamples,
                                       model=object@model,
                                       data=thisData)$value
                  
                  ## evaluate stopping rules
                  stopit <- stopTrial(object@stopping,
                                      dose=thisDose,
                                      samples=thisSamples,
                                      model=object@model,
                                      data=thisData)
                }
                
                # thisSeed<-get(".Random.seed")
                thisSeed<-simSeeds[iterSim]
                
                ## get the fit
                thisFit <- fit(object=thisSamples,
                               model=object@model,
                               data=thisData)
                
                ## DK: get the individuals stop reasons
                stopLength<-length(object@stopping@stopList)
                stop_sum <- matrix(nrow=1,ncol=stopLength)
                for (i in 1: stopLength){
                  stop_sum[,i]<- stopTrial(stopping=object@stopping@stopList[[i]], 
                                           dose=thisDose,
                                           samples=thisSamples, 
                                           model=object@model, 
                                           data=thisData)
                  
                }
                
                ## DK: get the median MTD and MAD
                thisMTDsamples<-dose(prob=object@nextBest@target,model=object@model,samples=thisSamples)
                thisMTDmedian <- median(thisMTDsamples)
                thisMTDMAD <- mad(thisMTDsamples)
                thisMTDCV<-(thisMTDMAD/thisMTDmedian)*100
                
                ## DK: get the true MTD from truth scenario if not provided by user
                if(is.null(trueMTD)){
                  targetDLT<- object@nextBest@target
                  thisTrueMTD <-
                    sapply(targetDLT,
                           function(t){
                             r <- try(uniroot(f=function(x){truth(x, ...) - t},
                                              lower=-1000000, upper=1000000)$root,
                                      silent=TRUE)
                             if(inherits(r, "try-error"))
                             {
                               return(NA_real_)
                             } else {
                               return(r)
                             }
                           })
                }else{
                  thisTrueMTD<- trueMTD
                }
                ## DK: calculate the median MTD relative error with regards to true MTD
                thisMTDmedianRE<-abs(1-thisMTDmedian/thisTrueMTD)*100
                
                
                ## return the results
                thisResult <-
                  list(data=thisData,
                       dose=thisDose,
                       fit=
                         subset(thisFit,
                                select=c(middle, lower, upper)),
                       stop=
                         attr(stopit,
                              "message"),
                       MTD=thisMTDmedian,
                       CV=thisMTDCV,
                       stop_sum=stop_sum,
                       trueMTD=thisTrueMTD,
                       MTDRE=thisMTDmedianRE,
                       stopNames=stopNames,
                       seed_iter=thisSeed
                       
                  )
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
                                              "mcmcOptions",
                                              "trueMTD"),
                                          parallel=if(parallel) nCores else NULL)
              
              ## put everything in the Simulations format:
              
              ## setup the list for the simulated data objects
              dataList <- lapply(resultList, "[[", "data")
              
              ## the vector of the final dose recommendations
              recommendedDoses <- as.numeric(sapply(resultList, "[[", "dose"))
              
              ## the vector of the final median MTD based on the model 
              medianMTD <- as.numeric(sapply(resultList, "[[", "MTD"))
              
              ## the vector of the final median CV% based on the model 
              medianCV <- as.numeric(sapply(resultList, "[[", "CV"))
              
              ## the vector of the final median MTD Rel. Error %  
              medianMTDRE <- as.numeric(sapply(resultList, "[[", "MTDRE"))
              
              ## the vector of the final true MTD based on the truth 
              trueScenarioMTD <- unique(as.numeric(sapply(resultList, "[[", "trueMTD")))
              
              
              ## setup the list for the final fits
              fitList <- lapply(resultList, "[[", "fit")
              
              ## the reasons for stopping
              stopReasons <- lapply(resultList, "[[", "stop")
              
              ## the reasons for stopping
              stopReasonsSum <- lapply(resultList, "[[", "stop_sum")
              
              ## the reasons for stopping
              seed_iterList <- as.numeric(sapply(resultList, "[[", "seed_iter"))
              
              
              
              ## return the results in the Simulations class object
              ret <- Simulations(data=dataList,
                                 doses=recommendedDoses,
                                 MTD=medianMTD,
                                 CV=medianCV,
                                 MTDRE=medianMTDRE,
                                 fit=fitList,
                                 stopReasons=stopReasons,
                                 stopReasonsSum=stopReasonsSum,
                                 trueMTD=trueScenarioMTD,
                                 seed=RNGstate,
                                 mcmcIterations=mcmcIterations,
                                 mcmcBurnin=mcmcBurnin,
                                 mcmcStep=mcmcStep,
                                 stopNames=stopNames,
                                 seed_iterSim=seed_iterList)
              
              return(ret)
            })
