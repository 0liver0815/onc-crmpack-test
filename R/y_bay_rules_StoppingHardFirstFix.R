#' Stopping rule based on the number of observed DLTs at a dose level.
#' This is a hards stopping rule. The unacceptable numer of DLTs are
#' specified in the protocol. I.e. in the protocol a certain numbers for 
#' stopping the study is specified
#' (often not based on a posterior distribution)
#'
#' @slot toxrule matrix with toxrule
#' @slot label label of stopping rule
#' 
#' @export
#' @keywords classes
.StoppingHardFirstFix <-
  setClass(Class="StoppingHardFirstFix",
           representation(toxrule="matrix",
                          label="character"),
           prototype(toxrule=matrix(c(3,2),nrow=2),
                     label="tox first hard fix"),
           contains="Stopping")


#' Initialization function for "StoppingHardFirstFix"
#'
#' @param toxrule see \code{\linkS4class{StoppingHardFirstFix}}
#' @param label see \code{\linkS4class{StoppingHardFirstFix}}
#'
#' @export
#' @keywords methods
StoppingHardFirstFix <- function(toxrule,
                                 label)
{
  .StoppingHardFirstFix(toxrule=toxrule,
                        label=label)
}

#' @describeIn stopTrial Stopping rule based on hard stopping rule
#' This method take into account the hardstopping rule, specified in the
#' protocol. I.e. in the protocol certain numbers for stopping are specified
#' (often not based on a posterior distribution)
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingHardFirstFix",
                      dose="numeric",
                      samples="Samples"),
          def=
            function(stopping, dose, samples, model, data, ...){
              # determine the doses that are considered to be toxic (exclude from escalation)
              dlttab <- table(data@y, data@x)
              if (sum(data@y)>0){
                firstdose <- colnames(dlttab)==data@doseGrid[1]
                if (sum(firstdose)>0){
                  #is the total number of patients at first dose matching to a rule?
                  matchN <- outer(stopping@toxrule[1,],sum(dlttab[,1]),"==")
                  #is the number of observed DLTs matching to a rule?
                  matchDLT <- outer(stopping@toxrule[2,],dlttab[rownames(dlttab)=="1",1],"<=")
                  #Is the total number and the number of DLTs matching
                  dosetox <- apply(matchN & matchDLT,2,sum)
                } else{dosetox <- 0} 
              } else{dosetox <- 0}
              
              # so can we stop?
              doStop <- dosetox > 0
              
              # generate message
              text <-
                paste(doStop, ':', stopping@label,
                      ": Hardstop rule that first dose is toxic",
                      ifelse(doStop, "do", "do not"),
                      "match, as",
                      dlttab[rownames(dlttab)=="1",1],
                      "DLTs occur in",
                      sum(dlttab[,1]),
                      "patients and thus",
                      ifelse(doStop, "equal or above", "below"),
                      "the maximum allowed number")
              
              # return both
              return(structure(doStop,
                               message=text))
            })
