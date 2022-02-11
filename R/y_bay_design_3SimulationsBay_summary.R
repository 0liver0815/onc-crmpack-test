######################################################################
# New class SimulationsBay to output some estimations

##' Class for the simulations output from model based designs
##'
##' This class captures the trial simulations from model based designs.
##' Additional slots fit and stopReasons compared to the general class
##' \code{\linkS4class{GeneralSimulations}}.
##'
##' @slot fit list with the final fits
##' @slot estimates list with estimates from the final run
##' @slot stopReasons list of stopping reasons for each simulation run
##'
##' @export
##' @keywords classes
.SimulationsBay <-
  setClass(Class="SimulationsBay",
           representation(fit="list",
                          stopReasons="list",
                          estimates="list"),
           ## note: this prototype is put together with the prototype
           ## for GeneralSimulations
           prototype(fit=
                       list(c(0.1, 0.2),
                            c(0.1, 0.2)),
                     stopReasons=
                       list("A", "A"),
                     estimates=
                       list(c(0.1, 0.2),
                            c(0.1, 0.2))),
           contains="GeneralSimulationsBay",
           validity=
             function(object){
               o <- Validate()
               
               nSims <- length(object@data)
               
               o$check(identical(length(object@fit), nSims),
                       "fit must have same length as data")
               o$check(identical(length(object@stopReasons), nSims),
                       "stopReasons must have same length as data")
               o$check(identical(length(object@estimates), nSims),
                       "estimates must have same length as data")
               
               o$result()
             })
validObject(.Simulations())


##' Initialization function for the "Simulations" class
##'
##' @param fit see \code{\linkS4class{Simulations}}
##' @param estimates see \code{\linkS4class{Simulations}}
##' @param stopReasons see \code{\linkS4class{Simulations}}
##' @param \dots additional parameters from \code{\link{GeneralSimulations}}
##' @return the \code{\linkS4class{Simulations}} object
##' @export
##'
##' @keywords methods
SimulationsBay <- function(fit,
                           estimates,
                           stopReasons,
                           ...)
{
  start <- GeneralSimulationsBay(...)
  .SimulationsBay(start,
                  fit=fit,
                  estimates=estimates,
                  stopReasons=stopReasons)
}



##' Summarize the model-based design simulations, relative to a given truth
##'
##' @param object the \code{\linkS4class{SimulationsBay}} object we want to
##' summarize
##' @param truth a function which takes as input a dose (vector) and returns the
##' true probability (vector) for toxicity
##' @param target the target toxicity interval (default: 30\%) used for the
##' computations
##' @param \dots Additional arguments can be supplied here for \code{truth}
##' @return an object of class \code{\linkS4class{SimulationsSummary}}
##'
##' @example examples/Simulations-method-summary.R
##' @export
##' @keywords methods
setMethod("summary",
          signature=
            signature(object="SimulationsBay"),
          def=
            function(object,
                     truth,
                     target=0.3,
                     ...){
              
              ## call the parent method
              start <- callNextMethod(object=object,
                                      truth=truth,
                                      target=target,
                                      ...)
              
              doseGrid <- object@data[[1]]@doseGrid
              
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
              
              ## extract MTD and CV 
              medianMTD <- sapply(object@estimates, function(l) l[[1]])
              cvMTD <-  sapply(object@estimates, function(l) l[[2]])
              
              ## Calculate Stop reasons in two ways
              
              ## extract label
              label <- strsplit(unlist(object@stopReasons[[1]]),' : ')
              stoplabel <- sapply(label, function(l) l[[2]])
              
              ## Count stop reasons in total
              ## split stop reason in two parts
              stop <- strsplit(unlist(object@stopReasons),' : ')
              stop <- lapply(stop, trimws)
              
              ## keep only first part and convert into logical vector
              stop1 <- as.logical(lapply(stop, function(l) l[[1]]))
              
              ## put vector into a matrix so that each row is equal to a study
              stop.res <- matrix(stop1,nrow=length(object@doses),byrow=T)
              
              ## calculate the stop reasons for each column of the matrix
              ## several stop reasons may occur at the same time
              ## the columns appear in the same order as they are put together in the simulation program
              ## myStoppinglow | myStoppinghigh | myStoppingCV | myStoppingnpat | myStoppingfirst
              stopreason <- apply(stop.res,2,sum)#/object@nsim*100
              names(stopreason) <- stoplabel
              
              ## count only the first occurrence, i.e. stop reason needs to be in
              ## sequential order of importance or logic
              stopreason2 <- apply(stop.res, 1, which)
              stopreason2 <- sapply(stopreason2, function(l) l[[1]])
              stopreason2 <- factor(stopreason2, levels=1:length(stoplabel))
              stopreason.unique <- table(stopreason2)#/object@nsim*100
              row.names(stopreason.unique) <- stoplabel              
              
              ## give back an object of class SimulationsSummary,
              ## for which we then define a print / plot method
              ret <- .SimulationsSummaryBay(
                start,
                fitAtDoseMostSelected=fitAtDoseMostSelected,
                meanFit=meanFit,
                medianMTD=medianMTD,
                cvMTD=cvMTD,
                stopreason=stopreason,
                stopreason.unique=stopreason.unique
              )
              
              return(ret)
            })


