#' Increments control based on the number of observed DLTs
#' This is a safety rule that limits further escalation based on the 
#' number of observed DLTs. The probability that a dose is toxic is calculated
#' using a Bin-Beta distribution with prior (a,b). If the probability exceeds
#' the threshold for a given dose, that dose and all doses above are excluded
#' from further escalation.
#' 
#' @slot target the target toxicity
#' @slot prob the threshold probability 
#' @slot prior.a the the first parameter of the prior Beta distribution
#' @slot prior.b the the second parameter of the prior Beta distribution
#' 
#' @export
#' @keywords classes 
.IncrementsSafetyStop <-
  setClass(Class="IncrementsSafetyStop",
           representation(target="numeric",
                          prob="numeric",
                          prior.a="numeric",
                          prior.b="numeric"),
           prototype(target=0.3,
                     prob=0.95,
                     prior.a=1,
                     prior.b=1),
           contains="Increments")

#' Initialization function for "IncrementsSafetyStop"
#'
#' @param target see \code{\linkS4class{IncrementsSafetyStop}}
#' @param prob see \code{\linkS4class{IncrementsSafetyStop}}
#' @param prior.a see \code{\linkS4class{IncrementsSafetyStop}}
#' @param prior.b see \code{\linkS4class{IncrementsSafetyStop}}
#'
#' @export
#' @keywords methods 
IncrementsSafetyStop <- function(target,
                                 prob,
                                 prior.a=1,
                                 prior.b=1)
{
  .IncrementsSafetyStop(target=target,
                        prob=prob,
                        prior.a=prior.a,
                        prior.b=prior.b)
}

#' @describeIn maxDose Determine the maximum possible dose for escalation
setMethod("maxDose",
          signature=
            signature(increments="IncrementsSafetyStop",
                      data="Data"),
          def=
            function(increments, data, ...){
              
              ## summary of observed data per dose level
              y <- factor(data@y, levels = c('0','1'))
              dlttab <- table(y, data@x)
              ## extract dose names as these get lost if only one dose available
              doses <- as.numeric(colnames(dlttab))
              
              ## toxicity probability per dose level
              toxprob <- 1-pbeta(increments@target, dlttab[2,] + increments@prior.a,
                                 apply(dlttab,2,sum) - dlttab[2,] + increments@prior.b)
              
              ## If return the min toxic dose level or maximum dose level if no dose is toxic
              if (sum(toxprob > increments@prob)>0){
                dosetox <- min(doses[which(toxprob > increments@prob)])
              }else{
                ## add 1 to max so that the max dose is always smaller
                dosetox <- max(data@doseGrid) +1
              } 
              
              ## determine doses that are not toxic
              doseok <- max(data@doseGrid[data@doseGrid < dosetox],0)
              
              return(doseok)
            })
