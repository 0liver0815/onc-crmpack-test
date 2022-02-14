
## -----------------------------------------------------------------------------------------
## Add pbeta in importFrom stats in the program crmPack-package.r
## -----------------------------------------------------------------------------------------

## --------------------------------------------------
## To be added in <[Rules-class.R]
## --------------------------------------------------

## -------------------------------------------------------
## Increments control based on Hard Safety Rule - Beta(a,b)
## -------------------------------------------------------

##' Increments control based on eligable dose levels according to Hard Safety
##' Rule - Beta(a,b)
##'
##' @slot target a probability > 0 and < 1 of the target toxicity rate
##' @slot prob a probability > 0 and < 1 of the confidence of the test
##' @slot a positive real number representing the shape parameter of a Beta(a,b) distribution
##' @slot b positive real number representing the shape parameter of a Beta(a,b) distribution
##'
##' @export
##' @keywords classes
.IncrementsHSRBeta_dk <-
  setClass(Class="IncrementsHSRBeta_dk",
           representation(target="numeric",
                          prob="numeric",
                          a="numeric",
                          b="numeric"),
           prototype(target=0.3,
                     prob=0.9,
                     a=1,
                     b=1),
           contains="Increments",
           validity=
             function(object){
               o <- Validate()

               o$check(is.probability(object@target,bounds=FALSE),
                       "target must be probability > 0 and < 1")
               o$check(is.probability(object@prob,bounds=FALSE),
                       "prob must be probability > 0 and < 1")
               o$check(is.numeric(object@a) & object@a > 0,
                       "Beta distribution shape parameter a must me a real number > 0")
               o$check(is.numeric(object@b) & object@b > 0,
                       "Beta distribution shape parameter b must me a real number > 0")

               o$result()
             })
validObject(.IncrementsHSRBeta_dk())

##' Initialization function for "IncrementsHSRBeta_dk"
##'
##' @param target see \code{\linkS4class{IncrementsHSRBeta_dk}}
##' @param prob see \code{\linkS4class{IncrementsHSRBeta_dk}}
##' @param a see \code{\linkS4class{IncrementsHSRBeta_dk}}
##' @param b see \code{\linkS4class{IncrementsHSRBeta_dk}}
##'
##' @return the \code{\linkS4class{IncrementsHSRBeta_dk}} object
##'
##' @export
##' @keywords methods
IncrementsHSRBeta_dk <- function(target,prob,a,b)
{
  .IncrementsHSRBeta_dk(target=target,
                     prob=prob,
                     a=a,
                     b=b)
}


## --------------------------------------------------
## To be added in <[Rules-methods.R]
## --------------------------------------------------

## --------------------------------------------------
## The maximum allowable number of dose levels method based on Hard Safety Rule - Beta(a,b)
## --------------------------------------------------

##' @describeIn maxDose Determine the maximum possible next dose based on
##' maximum eligible dose level to increment for the next dose
##' below the target toxicity, by a % probability based on the Beta(a,b) distribution.
##' this can be calculated for any number of patients in a cohort with a defined probability threshold of 90%
##' i.e. 1- pbeta(target, x+a, n-x+b) < 0.90
##' where x=DLTs, n=patients for a dose and use Beta(a=1,b=1) (or probably a=b=0.5)
##' if a dose satisfies the above criterio then eligible doses, are only the ones that that lower than that dose
##'
setMethod("maxDose",
          signature=
            signature(increments="IncrementsHSRBeta_dk",
                      data="Data"),
          def=
            function(increments, data, ...){
              #Load data dose vector and find the number of total subjects dosed at each dose
              n<-table(data@x)

              #Load data events vector at each dose and find the number of DLTs
              yd<-table(data@y,data@x)

              #check if both DLTs and no DLTs are present
              if (nrow(yd)==2){
                #Number of DLTs
                y1<-yd[2,]
              }

              #check if only DLTs or only no DLTs are present
              if (nrow(yd)==1){
                #Number of DLTs in each scenario:
                #if only DLTs
                if (dimnames(yd)[1]==1){y1=n}
                #if only no DLTs
                if (dimnames(yd)[1]==0){y1=0}
              }
              #Calculate the probability of DLTs at each dose  tested
              # i.e. 1- pbeta(target, y+a, n-y+b) < 0.90
              # where y=DLTs, n=patients for a dose and use Beta(a,b)
              probTarget<-(1- pbeta(increments@target, y1+increments@a, n-y1+increments@b))

              # find toxic doses
              toxic<-probTarget>=increments@prob

              # find safe doses
              safe<-probTarget<increments@prob

              # are there any safe doses
              safegrid<-sum(safe[]=="TRUE")

              # are there any toxic doses
              toxicgrid<-sum(toxic[]=="TRUE")

              # if only safe doses exist (toxic doses do not exist all) then all grid doses are eligible and max next dose is the maximum planned dose grid
              if (safegrid[]>0 && toxicgrid[]==0){
                maxdose<-max(data@doseGrid[])
              }else if (safegrid[]>0 && toxicgrid[]>0){
                # if safe and toxic doses both exist, then max next dose is the maximum safe dose
                maxdose<-data@doseGrid[safegrid]
              }else if (safegrid[]==0 && toxicgrid[]>0){
                # if safe doses do not exist and all doses toxic
                maxdose<-min(data@doseGrid)
                # no dose is eligible so we assign the lowest (first) dose. This will force the next recommended dose to be the lowest dose and will trigger the stopping rule stoppingruleHSFBeta.
                # therefore when using this increment rule it is necessary that you also use the relevant stoppingruleHSFBeta with the same setting. In this
                # way the absence of safe doses including the lowest one will trigger trial to stop as long as these 2 are set the same.

              }

              ## so the maximum next dose is
              ret <- maxdose

              return(ret)
            })

