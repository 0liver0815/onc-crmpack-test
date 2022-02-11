#' Increments control based on the number of observed DLTs
#' This is a safety rule that uses fixed numbers for observed DLTs
#' at a dose level to limit future escalation. I.e. In case that the number of DLTs
#' exceeds a predefined number specified in the protocol only lower doses or the
#' same dose are allowed for the further dose escalation.
#' 
#' Note: in case that samedose=TRUE, only the exact match of DLT is evaluated
#'       i.e. the function returns ok, if the number of observed DLTs are above
#'       the specified number. In this case the rule should only be used together
#'       with a rule using samedose=FALSE. The rule must be specified as a
#'       matrix with 2 rows.
#' 
#' @slot toxrule the toxrule for lower or the same dose. The rule must be given in a matrix form:
#' row1=Number of patients,
#' row2=number of DLTs
#' @slot samedose If true the same dose is ok
#' 
#' @export
#' @keywords classes 
.IncrementsSafetyStopFix <-
  setClass(Class="IncrementsSafetyStopFix",
           representation(toxrule="matrix",
                          samedose="logical"),
           prototype(toxrule=matrix(c(3,2,6,3),nrow=2),
                     samedose=FALSE),
           contains="Increments")



#' Initialization function for "IncrementsSafetyStopFix"
#'
#' @param toxrule see \code{\linkS4class{IncrementsSafetyStopFix}}
#' @param samedose see \code{\linkS4class{IncrementsSafetyStopFix}}
#'
#' @export
#' @keywords methods 
IncrementsSafetyStopFix <- function(toxrule,
                                    samedose=FALSE)
{
  .IncrementsSafetyStopFix(toxrule=toxrule,
                           samedose=samedose)
}

#' @describeIn maxDose Determine the maximum possible dose for escalation
setMethod("maxDose",
          signature=
            signature(increments="IncrementsSafetyStopFix",
                      data="Data"),
          def=
            function(increments, data, ...){
              # determine how many DLTs have occurred so far
              dltHappened <- sum(data@y)
              
              if (sum(data@y)>0){
                #create a table containing in the first row the number of patients without DLTs
                #and in the second row the number of patients with DLTs per dose
                dlttab <- table(data@y, data@x)
                
                #is the total number of patients at a dose matching to a value from toxrule?
                matchN <- outer(increments@toxrule[1,],apply(dlttab,2,sum),"==")
                
                #is the number of observed DLTs matching to a value from toxrule?
                if (increments@samedose){
                  matchDLT <- outer(increments@toxrule[2,],dlttab[rownames(dlttab)=="1",],"==")} # == instead of <=
                else{
                  matchDLT <- outer(increments@toxrule[2,],dlttab[rownames(dlttab)=="1",],"<=")}
                
                #Is the total number and the number of DLTs matching
                dosetox <- apply(matchN & matchDLT,2,sum)
                
                #If not, no dose must be excluded
                if (sum(dosetox)==0){doseok <- max(data@doseGrid)
                } else {
                  dosetox <- min(as.numeric(names(dosetox)[dosetox!=0]))
                  #exclude toxic dose and all doses above, or keep dose and exclude above
                  if (increments@samedose){
                    doseok  <- max(data@doseGrid[data@doseGrid <= dosetox],0)} # <= instead of <
                  else{
                    doseok  <- max(data@doseGrid[data@doseGrid < dosetox],0)}
                }
              }
              else{doseok <-max(data@doseGrid)}
              ret <- doseok
              
              return(ret)
            })


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

