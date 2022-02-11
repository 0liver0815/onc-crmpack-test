#' Hardstopping rule based on the Bin-Beta distribution and not on the
#' results from the model.The rule is triggered when
#' the first dose is considered to be toxic (above threshold prob) based on the
#' observed data and a Beta prior distribution. For calculation the posterior
#' bin-beta distribution, usually an uninformative Beta(1,1) prior is used.
#'
#' @slot target the target toxicity
#' @slot prob the threshold probability
#' @slot prior.a the first parameter of the prior Beta distribution
#' @slot prior.b the second parameter of the prior Beta distribution
#' @slot label label of stopping rule
#' 
#' @export
#' @keywords classes
.StoppingHardFirst <-
  setClass(Class="StoppingHardFirst",
           representation(target="numeric",
                          prob="numeric",
                          prior.a="numeric",
                          prior.b="numeric",
                          label="character"),
           prototype(target=0.33,
                     prob=0.9,
                     prior.a=1,
                     prior.b=1,
                     label="first tox hard dist"),
           contains="Stopping")

#' Initialization function for "StoppingHardFirst"
#'
#' @param target see \code{\linkS4class{StoppingHardFirst}}
#' @param prob see \code{\linkS4class{StoppingHardFirst}}
#' @param label see \code{\linkS4class{StoppingHardFirst}}
#' @param prior.a see \code{\linkS4class{StoppingHardFirst}}
#' @param prior.b see \code{\linkS4class{StoppingHardFirst}}
#' 
#' @export
#' @keywords methods
StoppingHardFirst <- function(target,
                              prob,
                              prior.a=1,
                              prior.b=1,
                              label)
{
  .StoppingHardFirst(target=target,
                     prob=prob,
                     prior.a=prior.a,
                     prior.b=prior.b,
                     label=label)
}


#' @describeIn stopTrial Stopping rule based on a hard stopping rule
#' This is a method that takes into account a hardstopping rule that is not
#' based on the model. The rule is triggered when the first dose is considered
#' to be toxic (above threshold prob) based on the observed data and an
#' uninformative Beta(1,1) prior (default).
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingHardFirst",
                      dose="numeric",
                      samples="Samples"),
          def=
            function(stopping, dose, samples, model, data, ...){
              # determine if the first doses is toxic
              # First dose Tested? 
              if (sum(data@x==data@doseGrid[1])>0){
                # summary data at first dose
                y <- factor(data@y, levels = c('0','1'))
                dlttab <- table(y, data@x)[,1]
                # Any DLT at first dose?
                if (dlttab[2]>0){
                  toxprobfirstdose <- 1-pbeta(stopping@target, dlttab[2]+stopping@prior.a,
                                              sum(dlttab)-dlttab[2]+stopping@prior.b)
                } else {toxprobfirstdose <- 0}
              } else {toxprobfirstdose <- 0}
              
              # so can we stop?
              doStop <- toxprobfirstdose > stopping@prob
              
              # generate message
              text <-
                paste(doStop, ' : ', stopping@label,
                      " : Hardstop probability of MTD below first dose is ",
                      round(toxprobfirstdose * 100),
                      "% and thus ",
                      ifelse(doStop, "above", "below"),
                      " the required ",
                      round(stopping@prob * 100),
                      "%",
                      sep="")
              
              # return both
              return(structure(doStop,
                               message=text))
            })
