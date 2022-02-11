#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Increments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


###########################################################################################
# A new increments class is introduced:
# Implementation of a safety rule that is based on bin-beta with uninformative beta(1,1) prior
# In case that the prob that a dose is toxic exceeds a predefined threshold, only lower doses
# are allowed in the further dose escalation
# The rule is implemented via the class Increments and named "IncrementsRelativeSafety"
# The related method is maxDose

##' Increments control based on relative safety
##'
##' @slot intervals a vector with the left bounds of the relevant intervals
##' @slot increments a vector of the same length with the maximum allowable
##' @slot target the target toxicity
##' @slot toxprob the threshold probability that the dose is above the target
##' 
##' @export
##' @keywords classes 
.IncrementsRelativeSafety <-
  setClass(Class="IncrementsRelativeSafety",
           representation(intervals="numeric",
                          increments="numeric",
                          target="numeric",
                          toxprob="numeric"),
           prototype(intervals=c(0, 2),
                     increments=c(2, 1),
                     traget=0.3,
                     toxprob=0.9),
           contains="Increments")



##' Initialization function for "IncrementsRelativeSafety"
##'
##' @param intervals see \code{\linkS4class{IncrementsRelativeSafety}}
##' @param increments see \code{\linkS4class{IncrementsRelativeSafety}}
##' @param target see \code{\linkS4class{IncrementsRelativeSafety}}
##' @param toxprob see \code{\linkS4class{IncrementsRelativeSafety}}
##'
##' @export
##' @keywords methods 
IncrementsRelativeSafety <- function(intervals,
                                     increments,
                                     target,
                                     toxprob)
{
  .IncrementsRelativeSafety(intervals=intervals,
                            increments=increments,
                            target=target,
                            toxprob=toxprob)
}


##' @describeIn maxDose Determine the maximum possible dose for escalation
##' @importFrom stats pbeta
setMethod("maxDose",
          signature=
            signature(increments="IncrementsRelativeSafety",
                      data="Data"),
          def=
            function(increments, data, ...){
              ## determine what was the last dose
              lastDose <- tail(data@x, 1)
              
              ## determine in which interval this dose was
              lastInterval <-
                findInterval(x=lastDose,
                             vec=increments@intervals)
              
              ## so the maximum next dose is
              ret <-
                (1 + increments@increments[lastInterval]) *
                lastDose
              
              ## determine the doses that are considered to be toxic (exclude from further dose escalation)
              ## the maximum allowed next non-toxic dose is stored in doseok
              if (sum(data@y)>0){
                dlttab <- table(data@y, data@x)
                toxprobdose <- 1-pbeta(increments@target, dlttab[rownames(dlttab)=="1",]+1, apply(dlttab,2,sum)-dlttab[rownames(dlttab)=="1",]+1)
                if (sum(toxprobdose>increments@toxprob)==0){doseok <- max(data@doseGrid)
                } else {
                  if (length(toxprobdose)==1){names(toxprobdose)<-colnames(dlttab)}
                  doseok <- max(data@doseGrid[data@doseGrid < min(as.numeric(names(toxprobdose)[toxprobdose>increments@toxprob]))],0)}
              }else {doseok<-max(data@doseGrid)}
              
              ##return the mimum dose of the increments and the safe doses
              ret <- min(ret,doseok)
              
              return(ret)
            })


###########################################################################################
# Implementation of a safety rule that uses fixed numbers from the protocol
# In case that the number of DLTs exceeds a predefined number from the protocol
# only lower or the same dose are allowed in the further dose escalation
# only lower dose: toxrule
# only lower and the same: toxrulesame
# The rule is implemented via the class Increments and named "IncrementsRelativeSafetyFix"
# The rule must be specified as a matrix: row1=Number of patients, row2=number of tox
# The related method is maxDose

##' Increments control based on relative safety
##'
##' @slot intervals a vector with the left bounds of the relevant intervals
##' @slot increments a vector of the same length with the maximum allowable
##' @slot toxrule the toxrule for lower doses only
##' @slot toxrulesame toxrule for lower and the same dose dose
##' 
##' @export
##' @keywords classes 
.IncrementsRelativeSafetyFix <-
  setClass(Class="IncrementsRelativeSafetyFix",
           representation(intervals="numeric",
                          increments="numeric",
                          toxrule="matrix",
                          toxrulesame="matrix"),
           prototype(intervals=c(0, 2),
                     increments=c(2, 1),
                     toxrule=matrix(c(3,2,6,3),nrow=2),
                     toxrulesame=matrix(c(3,1,6,2),nrow=2)),
           contains="Increments")



##' Initialization function for "IncrementsRelativeSafetyFix"
##'
##' @param intervals see \code{\linkS4class{IncrementsRelativeSafetyFix}}
##' @param increments see \code{\linkS4class{IncrementsRelativeSafetyFix}}
##' @param toxrule see \code{\linkS4class{IncrementsRelativeSafetyFix}}
##' @param toxrulesame see \code{\linkS4class{IncrementsRelativeSafetyFix}}
##'
##' @export
##' @keywords methods 
IncrementsRelativeSafetyFix <- function(intervals,
                                        increments,
                                        toxrule,
                                        toxrulesame)
{
  .IncrementsRelativeSafetyFix(intervals=intervals,
                               increments=increments,
                               toxrule=toxrule,
                               toxrulesame=toxrulesame)
}

##' @describeIn maxDose Determine the maximum possible dose for escalation
setMethod("maxDose",
          signature=
            signature(increments="IncrementsRelativeSafetyFix",
                      data="Data"),
          def=
            function(increments, data, ...){
              ## determine what was the last dose
              #lastDose <- tail(data@x, 1)
              # determine what was the maximum tested dose
              maxDose <- max(data@x)
              
              ## determine in which interval this dose was
              maxInterval <-
                findInterval(x=maxDose,
                             vec=increments@intervals)
              
              ## so the maximum next dose is
              ret <-
                (1 + increments@increments[maxInterval]) * maxDose
              
              ## determine the doses that are considered to be too toxic (exclude from further dose escalation)
              ## the maximum allowed next non-toxic dose is stored in doseok
              if (sum(data@y)>0){
                dlttab <- table(data@y, data@x)
                #is the total number of patients at a dose matching to a rule?
                matchN <- outer(increments@toxrule[1,],apply(dlttab,2,sum),"==")
                #is the number of observed DLTs matching to a rule?
                matchDLT <- outer(increments@toxrule[2,],dlttab[rownames(dlttab)=="1",],"<=")
                #Is the total number and the number of DLTs matching
                dosetox <- apply(matchN & matchDLT,2,sum)
                #If not, no dose must be excluded
                if (sum(dosetox)==0){doseok <- max(data@doseGrid)
                } else {
                  #exclude all doses above the lowest toxic dose
                  dosetox <- min(as.numeric(names(dosetox)[dosetox!=0]))
                  doseok <- max(data@doseGrid[data@doseGrid < dosetox],0)
                }
                
                if (sum(increments@toxrulesame) != 0){
                  matchNsame <- outer(increments@toxrulesame[1,],apply(dlttab,2,sum),"==")
                  matchDLTsame <- outer(increments@toxrulesame[2,],dlttab[rownames(dlttab)=="1",],"==") # == instead of <=
                  dosetoxsame <- apply(matchNsame & matchDLTsame,2,sum)
                  if (sum(dosetoxsame)==0){doseoksame <- max(data@doseGrid)
                  } else {
                    #exclude all doses above the lowest toxic dose (include the dose with DLTs specified in toxrulesame)
                    dosetoxsame <- min(as.numeric(names(dosetoxsame)[dosetoxsame!=0]))
                    doseoksame <- max(data@doseGrid[data@doseGrid <= dosetoxsame],0) # <= instead of <
                  }
                }
                else {doseoksame <-max(data@doseGrid)}
                
              } else {doseok <- doseoksame <- max(data@doseGrid)}
              
              
              ##return the maximum dose allowed in combination with the safety rules
              ret <- min(ret, doseok, doseoksame)
              
              return(ret)
            })


###########################################################################################
# Use the existing IncrementsNumDoseLevels and adds rule for fixed numbers of DLTs:
# Add the implementation of a safety rule that uses fixed numbers from the protocol
# In case that the number of DLTs exceeds a predefined number from the protocol
# only lower or the same dose are allowed in the further dose escalation
# only lower dose: toxrule
# only lower and the same: toxrulesame
# The rule is implemented via the class Increments and named "IncrementsNumDoseLevelsBaySafetyFix"
# The rule must be specified as a matrix: row1=Number of patients, row2=number of toxicities
# e.g. toxrule=matrix(c(3,3,4,3,6,4,7,5,9,6,10,6),nrow=2) 3/3, 4/3 etc.
# The related method is maxDose
#######################################################################################################

##' Increments control based the maximum dose level increase
##' (not relative to the dose, absolute to the number of the dose level)
##' and a rule spefied in the protocol
##'
##' @slot maxLevels a vector with the left bounds of the relevant intervals
##' @slot toxrule the toxrule for lower doses only
##' @slot toxrulesame toxrule for lower and the same dose dose
##' 
##' @export
##' @keywords classes
.IncrementsNumDoseLevelsBaySafetyFix <-
  setClass(Class="IncrementsNumDoseLevelsBaySafetyFix",
           representation(maxLevels="integer",
                          toxrule="matrix",
                          toxrulesame="matrix"),
           prototype(maxLevels=1L,
                     toxrule=matrix(c(3,2,6,3),nrow=2),
                     toxrulesame=matrix(c(3,1,6,2),nrow=2)),
           contains="Increments")

##' Initialization function for "IncrementsRelativeSafetyFix"
##'
##' @param maxLevels see \code{\linkS4class{IncrementsNumDoseLevelsBaySafetyFix}}
##' @param toxrule see \code{\linkS4class{IncrementsNumDoseLevelsBaySafetyFix}}
##' @param toxrulesame see \code{\linkS4class{IncrementsNumDoseLevelsBaySafetyFix}}
##'
##' @export
##' @keywords methods 
IncrementsNumDoseLevelsBaySafetyFix <- function(maxLevels=1,
                                                toxrule,
                                                toxrulesame)
{
  .IncrementsNumDoseLevelsBaySafetyFix(maxLevels=safeInteger(maxLevels),
                                       toxrule=toxrule,
                                       toxrulesame=toxrulesame)
}

##' @describeIn maxDose Determine the maximum possible dose for escalation
setMethod("maxDose",
          signature=
            signature(increments="IncrementsNumDoseLevelsBaySafetyFix",
                      data="Data"),
          def=
            function(increments, data, ...){
              ## determine what was the level of the last dose
              maxDoseLevel <- max(data@xLevel)
              
              ## determine the maximum next dose level
              maxNextDoseLevel <- min(length(data@doseGrid),
                                      maxDoseLevel + increments@maxLevels)
              
              ## so the maximum next dose is
              ret <- data@doseGrid[maxNextDoseLevel]
              
              ## determine the doses that are considered to be too toxic (exclude from further dose escalation)
              ## the maximum allowed next non-toxic dose is stored in doseok
              if (sum(data@y)>0){
                dlttab <- table(data@y, data@x)
                #is the total number of patients at a dose matching to a rule?
                matchN <- outer(increments@toxrule[1,],apply(dlttab,2,sum),"==")
                #is the number of observed DLTs matching to a rule?
                matchDLT <- outer(increments@toxrule[2,],dlttab[rownames(dlttab)=="1",],"<=")
                #Is the total number and the number of DLTs matching
                dosetox <- apply(matchN & matchDLT,2,sum)
                #If not, no dose must be excluded
                if (sum(dosetox)==0){doseok <- max(data@doseGrid)
                } else {
                  #exclude all doses above the lowest toxic dose
                  dosetox <- min(as.numeric(names(dosetox)[dosetox!=0]))
                  doseok <- max(data@doseGrid[data@doseGrid < dosetox],0)
                }
                
                if (sum(increments@toxrulesame) != 0){
                  matchNsame <- outer(increments@toxrulesame[1,],apply(dlttab,2,sum),"==")
                  matchDLTsame <- outer(increments@toxrulesame[2,],dlttab[rownames(dlttab)=="1",],"==") # == instead of <=
                  dosetoxsame <- apply(matchNsame & matchDLTsame,2,sum)
                  if (sum(dosetoxsame)==0){doseoksame <- max(data@doseGrid)
                  } else {
                    #exclude all doses above the lowest toxic dose (include the dose with DLTs specified in toxrulesame)
                    dosetoxsame <- min(as.numeric(names(dosetoxsame)[dosetoxsame!=0]))
                    doseoksame <- max(data@doseGrid[data@doseGrid <= dosetoxsame],0) # <= instead of <
                  }
                }
                else {doseoksame <-max(data@doseGrid)}
                
              } else {doseok <- doseoksame <- max(data@doseGrid)}
              
              ##return the minimum dose of the increments and the safe doses
              ret <- min(ret, doseok, doseoksame)
              
              return(ret)
            })


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# Obsolete Stopping rules
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#Change existing stopping rule:
## --------------------------------------------------
## Stopping based on number of patients near to next best dose
## --------------------------------------------------

##' Stop based on number of patients near to next best dose
##'
##' @slot nPatients number of required patients
##' @slot percentage percentage (between 0 and 100) within the next best dose
##' the patients must lie
##' 
##' @example examples/Rules-class-StoppingPatientsNearDose.R
##' @keywords classes
##' @export
.StoppingPatientsNearDoseBay <-
  setClass(Class="StoppingPatientsNearDoseBay",
           representation(nPatients="integer",
                          percentage="numeric"),
           prototype(nPatients=10L,
                     percentage=50),
           contains="Stopping",
           validity=function(object){
             o <- Validate()
             
             o$check((object@nPatients > 0L) && is.scalar(object@nPatients),
                     "nPatients must be positive scalar")
             o$check(is.probability(object@percentage / 100),
                     "percentage must be between 0 and 100")
             
             o$result()
           })
validObject(.StoppingPatientsNearDoseBay())


##' Initialization function for "StoppingPatientsNearDose"
##'
##' @param nPatients see \code{\linkS4class{StoppingPatientsNearDose}}
##' @param percentage see \code{\linkS4class{StoppingPatientsNearDose}}
##' @return the \code{\linkS4class{StoppingPatientsNearDose}} object
##'
##' @export
##' @keywords methods
StoppingPatientsNearDoseBay <- function(nPatients,
                                        percentage)
{
  .StoppingPatientsNearDoseBay(nPatients=safeInteger(nPatients),
                               percentage=percentage)
} 

## -------------------------------------------------------------
## Stopping based on number of patients near to next best dose
## -------------------------------------------------------------

##' @describeIn stopTrial Stop based on number of patients near to next best
##' dose
##' 
##' @example examples/Rules-method-stopTrial-StoppingPatientsNearDose.R
setMethod("stopTrial",
          signature=
            signature(stopping="StoppingPatientsNearDoseBay",
                      dose="numeric",
                      samples="ANY",
                      model="ANY",
                      data="Data"),
          def=
            function(stopping, dose, samples, model, data, ...){
              ## determine the range where the cohorts must lie in
              lower <- (100 - stopping@percentage) / 100 * dose
              upper <- (100 + stopping@percentage) / 100 * dose
              
              ## how many patients lie there?
              nPatients <- sum((data@x >= lower) & (data@x <= upper))
              
              ## so can we stop?
              doStop <- nPatients >= stopping@nPatients
              
              ## generate message
              text <- paste(doStop, ' : ', nPatients,
                            " patients lie within ",
                            stopping@percentage,
                            "% of the next best dose ",
                            dose,
                            ". This ",
                            ifelse(doStop, "reached", "is below"),
                            " the required ",
                            stopping@nPatients,
                            " patients",
                            sep="")
              
              ## return both
              return(structure(doStop,
                               message=text))
            }) 


## -------------------------------------------------------------
## Next best dose specific for a model
## -------------------------------------------------------------


##' The class with the input for finding the next best MTD estimate with highest
##' probability to be the MTD
##'
##' @slot target the target toxicity probability
##' 
##' @export
##' @keywords classes
.NextBestKadaneBay <- setClass(Class = "NextBestKadaneBay",
                               contains = "NextBest", representation(target = "numeric"))


##' Initialization function for class "NextBestKadaneBay"
##'
##' @param target see \code{\linkS4class{NextBestKadaneBay}}
##'
##' @export
##' @keywords methods
NextBestKadaneBay <- function(target) {
  .NextBestKadaneBay(target = target)
}

##' @describeIn nextBest Find the next best dose based on MTD estimate with highest
##' probability to be the MTD
setMethod("nextBest",
          signature = signature(nextBest = "NextBestKadaneBay",
                                doselimit = "numeric", samples = "Samples", model = "Model",
                                data = "Data"),
          def = function(nextBest, doselimit, samples, model, data, ...) {
            #determine maximum allowed index
            maxdoseOK <- if (length(doselimit)) {
              max(which(data@doseGrid <= doselimit),0)
            } else {data@nGrid}
            # create matrix that repeats the dosegrid per row
            probDose <- matrix(data@doseGrid, nrow=length(samples@data$p0), ncol=data@nGrid, byrow=T)
            # calculate the toxicity per dose per sample
            probdose <- 1/(1 + exp(-logit(samples@data$p0) - ((logit(nextBest@target)-logit(samples@data$p0))/samples@data$MTD) * probDose))
            #count the number of dose levels below target=next best dose for iteration
            #and calculate the frequencies of the recommendation
            probMTDdist <- table(apply(probdose<nextBest@target,1,sum))
            #determine the dose with the highest frequency
            bestIndex <- as.integer(names(which.max(probMTDdist)))
            bestIndex <- min(bestIndex, maxdoseOK)
            if (bestIndex==0){bestDose<-data@doseGrid[1]
            }else{bestDose <- data@doseGrid[bestIndex]}
            return(list(value = bestDose))
          })


##' The class with the input for finding the next best MTD estimate with highest
##' probability to be the MTD
##'
##' @slot target the target toxicity probability
##' 
##' @export
##' @keywords classes
.NextBestTwoParBay <- setClass(Class = "NextBestTwoParBay",
                               contains = "NextBest",
                               representation(target = "numeric"))


##' Initialization function for class "NextBestTwoParBay"
##'
##' @param target see \code{\linkS4class{NextBestTwoParBay}}
##'
##' @export
##' @keywords methods
NextBestTwoParBay <- function(target) {
  .NextBestTwoParBay(target = target)
}

##' @describeIn nextBest Find the next best dose based on MTD estimate with highest
##' probability to be the MTD
setMethod("nextBest",
          signature = signature(nextBest = "NextBestTwoParBay",
                                doselimit = "numeric",
                                samples = "Samples",
                                model = "Model",
                                data = "Data"),
          def = function(nextBest, doselimit, samples, model, data, ...) {
            
            #determine maximum allowed index
            maxdoseOK <- if (length(doselimit)) {
              max(which(data@doseGrid <= doselimit),0)
            } else {data@nGrid}
            
            # create empty matrix
            probDose <- matrix(data@doseGrid, nrow=length(samples@data$int), ncol=data@nGrid, byrow=T)
            
            # calculate the toxicity per dose per sample
            probdose <- 1/(1+exp(-samples@data$int-samples@data$slope * probDose))
            
            #count the number of dose levels below target=next best dose for iteration
            #and calculate the frequencies of the recommendation
            probMTDdist <- table(apply(probdose<nextBest@target,1,sum))
            
            #determine the dose with the highest frequency
            bestIndex <- as.integer(names(which.max(probMTDdist)))
            bestIndex <- min(bestIndex, maxdoseOK)
            if (bestIndex==0){bestDose<-data@doseGrid[1]
            }else{bestDose <- data@doseGrid[bestIndex]}
            
            ## determine what was the last dose
            lastDose <- tail(data@x, 1)
            
            ## determine number of patients on last dose and number of DLTs
            lastDose.n <- length(data@y[data@cohort==tail(data@cohort,1)])
            lastDose.y <- sum(data@y[data@cohort==tail(data@cohort,1)])
            
            ## if a DLT occurs and the cohort size is equal to 1 expand at the same dose level
            ## regardless of any other result from increments, i.e. overwrite individual results
            if (lastDose.y > 0 & lastDose.n == 1){bestDose <- lastDose}
            
            return(list(value = bestDose,
                        allocation=probMTDdist))
          })

