#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
#
# Increments
#
#
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## --------------------------------------------------
## Increments control based on the number of patients at the last dose level and DLTs
## This is a special case for the ATRi protocol and introduced to mitigate the problem
## that if 1 DLT in 1 subject occur, the cohort should be expanded to 3 subjects. This 
## is not easily possible with the current mechanism in crmpack
## --------------------------------------------------
## This controls the increment in case of small cohort size, so that no increase based on the model
## is possible, i.e. if cohort size is 1 and 1 patient is observed with one DLT
## the dose should only be the same as, i.e. the cohort should be expanded to the full
## cohort size

##' Increments control based on relative differences in terms of DLTs
##'
##' Note that \code{DLTintervals} is to be read as follows. If for example,
##' we want to specify three intervals: First 0 DLTs, second 1 or 2 DLTs, and
##' third at least 3 DLTs, then we specify
##' \code{DLTintervals} to be \code{c(0, 1, 3)}. That means, the right
##' bound of the intervals are exclusive to the interval -- the vector only
##' gives the left bounds of the intervals. The last interval goes from 3 to
##' infinity.
##'
##' @slot NlastDL an integer with the min number of patients in case that a DLT occurs
##' @slot DLTlastDL an integer with the number of DLTs where an expansion should happen  
##' relative increments in the \code{DLTintervals}
##'
## @example examples/Rules-class-IncrementsNlastDL.R
##' @export
##' @keywords classes
.IncrementsNlastDL <-
  setClass(Class="IncrementsNlastDL",
           representation(NlastDL="integer",
                          DLTlastDL="integer"),
           prototype(NlastDL=as.integer(2),
                     DLTlastDL=as.integer(1)),
           contains="Increments")

##' Initialization function for "IncrementsRelativeDLT"
##'
##' @param NlastDL see \code{\linkS4class{IncrementsNlastDL}}
##' @param DLTlastDL see \code{\linkS4class{IncrementsNlastDL}}
##' @return the \code{\linkS4class{IncrementsNlastDL}} object
##'
##' @export
##' @keywords methods
IncrementsNlastDL <- function(NlastDL,
                              DLTlastDL)
{
  .IncrementsNlastDL(NlastDL=safeInteger(NlastDL),
                     DLTlastDL=safeInteger(DLTlastDL))
}

## --------------------------------------------------
## The maximum allowable relative increments in terms of DLTs
## --------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' relative increments determined by DLTs so far
##' 
## @example examples/Rules-method-maxDose-IncrementsNlastDL.R
setMethod("maxDose",
          signature=
            signature(increments="IncrementsNlastDL",
                      data="Data"),
          def=
            function(increments, data, ...){
              ## determine what was the last dose
              lastDose <- tail(data@x, 1)
              
              ## determine number of patients on last dose and number of DLTs
              lastDose.n <- length(data@y[data@cohort==tail(data@cohort,1)])
              lastDose.y <- sum(data@y[data@cohort==tail(data@cohort,1)])
              
              ## if at least the number of DLTs occur and the cohort size is below the specified number
              ## no escalation should be possible
              if (lastDose.y >= increments@DLTlastDL & lastDose.n < increments@NlastDL){ret <- lastDose}
              else{ret <- data@doseGrid[data@nGrid]}

              return(ret)
            })


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## -----------------------------------------------------------
## The maximum allowable relative increments taking the first observed DLT into account
## -----------------------------------------------------------

##' Max increment based on minimum of multiple increment rules
##'
##' This class can be used to combine multiple increment rules with the MIN
##' operation.
##'
##' \code{IncrementsList} contains all increment rules, which are again
##' objects of class \code{\linkS4class{Increments}}. The minimum of these
##' individual increments is taken to give the final maximum increment.
##'
##' @slot IncrementsList list of increment rules
##'
##' @example examples/Rules-class-IncrementMin.R
##' @keywords classes
##' @export
.IncrementMinDLT1N1 <-
  setClass(Class="IncrementMinDLT1N1",
           representation(IncrementsList="list"),
           prototype(IncrementsList=
                       list(IncrementsRelativeDLT(DLTintervals=as.integer(c(0, 1)),
                                                  increments=c(2, 1)),
                            IncrementsRelative(intervals=c(0, 2),
                                               increments=c(2, 1)))),
           contains="Increments")


##' Initialization function for "IncrementMinDLT1N1"
##'
##' @param IncrementsList see \code{\linkS4class{IncrementMinDLT1N1}}
##' @return the \code{\linkS4class{IncrementMinDLT1N1}} object
##'
##' @export
##' @keywords methods
IncrementMinDLT1N1 <- function(IncrementsList)
{
  .IncrementMinDLT1N1(IncrementsList=IncrementsList)
}


## --------------------------------------------------
## The maximum allowable relative increments taking the first observed DLT into account
## --------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' multiple increment rules (taking the minimum across individual increments).
##' 
## @example examples/Rules-method-maxDose-IncrementMinDLT1N1.R
setMethod("maxDose",
          signature=
            signature(increments="IncrementMinDLT1N1",
                      data="Data"),
          def=
            function(increments, data, ...){
              
              ## apply the multiple increment rules
              individualResults <-
                sapply(increments@IncrementsList,
                       maxDose,
                       data=data,
                       ...)

              ## so the maximum increment is the minimum across the individual increments 
              ret <- min(individualResults)
              
              ## determine what was the last dose
              lastDose <- tail(data@x, 1)
              
              ## determine number of patients on last dose and number of DLTs
              lastDose.n <- length(data@y[data@cohort==tail(data@cohort,1)])
              lastDose.y <- sum(data@y[data@cohort==tail(data@cohort,1)])
              
              ## if a DLT occurs and the cohort size is equal to 1 expand at the same dose level
              ## regardless of any other result from increments, i.e. overwrite individual results
              if (lastDose.y > 0 & lastDose.n == 1){ret <- lastDose}

              return(ret)
            })

