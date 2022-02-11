
## -----------------------------------------------------------------------------------------
## Add geom_rect in importFrom ggplot2 in the program crmPack-package.r
## -----------------------------------------------------------------------------------------

## -----------------------------------------------------------------------------------------
## New addition for [Rules-class.R] 
## -----------------------------------------------------------------------------------------

## ------------------------------------------------------
## Next best dose based on CRM MTD allocation probability
## ------------------------------------------------------

##' The class with the input for finding the next best dose based on CRM MTD 
##' allocation probability.
##'
##' @slot target the target toxicity probability
##' @slot pbomethod the method used to handle the allocation probability for the
##' zero dose. Applicable to cases where data@plabebo = "FALSE"  
##' If \code{none}, then the number of simulation instances allocated to zero 
##' dose, are exuded from the derivation of allocation probabilities.
##' If \code{min}, then the number of simulation instances allocated to zero 
##' dose, are combined with the instances allocated to the minimum planned dose.
##' If \code{max}, then the number of simulation instances allocated to zero 
##' dose, are combined with the instances allocated to the maximum planned dose.
##' 
##' @export
##' @keywords classes
.NextBestMTDCRM <-
  setClass(Class="NextBestMTDCRM",
           representation(target="numeric",
                          pbomethod="character"),
           prototype(target=0.3,
                     pbomethod="none"),
           contains=list("NextBest"),
           validity=
             function(object){
               o <- Validate()
               
               o$check(is.probability(object@target,
                                      bounds=FALSE),
                       "target must be probability > 0 and < 1")
               
               o$check(is.scalar(object@pbomethod) && object@pbomethod %in% c("none", "min", "max"),
                       "pbomethod must be either 'none' or 'min' or 'max'")
               
               o$result()
             })
validObject(.NextBestMTDCRM())

##' Initialization function for class "NextBestMTDCRM"
##'
##' @param target see \code{\linkS4class{NextBestMTDCRM}}
##' @param pbomethod see \code{\linkS4class{NextBestMTDCRM}}
##' @return the \code{\linkS4class{NextBestMTDCRM}} object
##'
##' @export
##' @keywords methods

NextBestMTDCRM <- function(target,
                            pbomethod=c("none", "min" , "max"))
{
  pbomethod <- match.arg(pbomethod)
  .NextBestMTDCRM(target=target,
                   pbomethod=pbomethod)
}

## -----------------------------------------------------------------------------------------
## New addition for [Rules-methods.R] 
## -----------------------------------------------------------------------------------------


## --------------------------------------------------
## The CRM MTD method
## --------------------------------------------------

##' @describeIn nextBest Find the next best dose based on CRM MTD allocation 
##' probability
##'
##' @importFrom ggplot2 ggplot geom_density xlab ylab xlim aes geom_vline
##' geom_text

setMethod("nextBest",
          signature=
            signature(nextBest="NextBestMTDCRM",
                      doselimit="numeric",
                      samples="Samples",
                      model="Model",
                      data="Data"),
          def=
            function(nextBest, doselimit, samples, model, data, ...){
              
              if(identical(length(doselimit), 0L))
              {
                warning("doselimit is empty, therefore no dose limit will be applied")
              }
              
              ## first we have to get samples from the dose-tox
              ## curve at the dose grid points.
              probSamples <- matrix(nrow=sampleSize(samples@options),
                                    ncol=data@nGrid)
              
              ## evaluate the probs, for all samples.
              for(i in seq_len(data@nGrid))
              {
                ## Now we want to evaluate for the
                ## following dose:
                probSamples[, i] <- prob(dose=data@doseGrid[i],
                                         model,
                                         samples)
              }
              
              ## count the total number of samples 
              nrow <- sampleSize(samples@options)
              
              ## Now identify dose rate probabilities that are less than or 
              ## equal to target 
              probTarget <-(probSamples <= nextBest@target)
              
              ## sum the number of doses for each row, that meet the 
              ## criterion (TRUE) 
              sumMTD <- rowSums(probTarget[]=='TRUE')
              
              ## count the number of different values of the sums 
              countsumMTD <- table(sumMTD)
              
              ## if the dimension of  countsumMTD ias not the same as the nGrid+1
              ## we need to impute with zero's the missing dose levels
              if (dim(countsumMTD)<(data@nGrid+1)){
                frame<-as.data.frame(countsumMTD)
                framefull<-data.frame(sumMTD=0:data@nGrid)
                frame1 <- merge(x=framefull ,y= frame,by="sumMTD", all.x = TRUE)
                frame1$Freq[is.na(frame1$Freq)] = 0
                countsumMTD <- as.table(frame1$Freq)
                names(countsumMTD)<-c(0:data@nGrid)
              }
              
              ## percentage of different values of the sums 
              freqMTD <- (countsumMTD/nrow)*100
              
              ## when the data object definition does NOT contain a placebo dose
              if(data@placebo == "FALSE")
              {
                
                ## if raw frequencies corresponding to planned doses are used for the nexbest dose
                if(nextBest@pbomethod == "none")
                {
                  countsumMTDnone<-countsumMTD[2:dim(countsumMTD)]
                  rr <- (countsumMTDnone/sum(countsumMTDnone))*100
                  rrmaxvalue<-rr[rr[]==max(rr)]
                  if (length(rrmaxvalue)==1){
                    doselevelMTD<-min(as.numeric(rownames(data.frame(rrmaxvalue))))
                    doseMTD<-data@doseGrid[doselevelMTD]
                  }
                  ## in case of ties in rrmaxvalue's the lowest index/dose is selected 
                  if (length(rrmaxvalue)>1){
                    doselevelMTD<-min(as.numeric(rownames(rrmaxvalue)))
                    doseMTD<-data@doseGrid[doselevelMTD]
                  }
                }
                
                ## if cumulative frequencies corresponding to planned doses are used for the nexbest dose
                ## and the frequency of zero dose is added to the frequency of the lowest planned dose
                if(nextBest@pbomethod == "min")
                {
                  freqMTDdmin <-freqMTD
                  freqMTDdmin[2]<-freqMTDdmin[1]+freqMTDdmin[2]
                  rr<-freqMTDdmin[2:dim(freqMTDdmin)]
                  rrmaxvalue<-rr[rr[]==max(rr)]
                  if (length(rrmaxvalue)==1){
                    doselevelMTD<-min(as.numeric(rownames(data.frame(rrmaxvalue))))
                    doseMTD<-data@doseGrid[doselevelMTD]
                  }
                  ## in case of ties in rrmaxvalue's the lowest index/dose is selected 
                  if (length(rrmaxvalue)>1){
                    doselevelMTD<-min(as.numeric(rownames(rrmaxvalue)))
                    doseMTD<-data@doseGrid[doselevelMTD]
                  }
                }
                
                ## if cumulative frequencies corresponding to planned doses are used for the nexbest dose
                ## and the frequency of zero dose is added to the frequency of the highest planned dose
                if(nextBest@pbomethod == "max")
                {
                  freqMTDdmax <-freqMTD
                  freqMTDdmax[dim(freqMTDdmax)]<-freqMTDdmax[1]+freqMTDdmax[dim(freqMTDdmax)]
                  rr<-freqMTDdmax[2:dim(freqMTDdmax)] 
                  rrmaxvalue<-rr[rr[]==max(rr)]
                }
                if (length(rrmaxvalue)==1){
                  doselevelMTD<-min(as.numeric(rownames(data.frame(rrmaxvalue))))
                  doseMTD<-data@doseGrid[doselevelMTD]
                }
                ## in case of ties in rrmaxvalue's the lowest index/dose is selected 
                if (length(rrmaxvalue)>1){
                doselevelMTD<-min(as.numeric(rownames(rrmaxvalue)))
                doseMTD<-data@doseGrid[doselevelMTD]
                }
                
              }
              
              ## when the data object definition contains a placebo dose
              ## if TRUE the first dose level in the grid is considered as PLACEBO
              if(data@placebo == "TRUE")
              {
                rr<-freqMTD
                rrmaxvalue<-rr[rr[]==max(rr)]
                if (length(rrmaxvalue)==1){
                  doselevelMTD<-min(as.numeric(rownames(data.frame(rrmaxvalue))))
                  doseMTD<-data@doseGrid[doselevelMTD]
                }
                ## in case of ties in rrmaxvalue's the lowest index/dose is selected 
                if (length(rrmaxvalue)>1){
                  doselevelMTD<-min(as.numeric(rownames(rrmaxvalue)))
                  doseMTD<-data@doseGrid[doselevelMTD]
                }
              }
              
              ## be sure which doses are ok with respect to maximum
              ## possible dose - if one was specified
              dosesOK <-
                if(length(doselimit)){
                  which(data@doseGrid <= doselimit)
                }else{
                seq_along(data@doseGrid)
                }
              
              dosesOKmax<-max(data@doseGrid[dosesOK])
              
              if (doseMTD >= dosesOKmax){
                ret <- dosesOKmax
              }
              
              if (doseMTD < dosesOKmax) {
                ret <- doseMTD
              }

              ## produce plot
              plot1 <- ggplot() +
                geom_line(data=data.frame(x=data@doseGrid,
                                          y=as.numeric(rr) ),
                          aes(x=x, y=y),
                          linetype = "dashed",
                          size=1,
                          colour="red"
                ) +
                geom_point(data=data.frame(x=data@doseGrid,
                                           y=as.numeric(rr) 
                ),
                aes(x=x, y=y),
                size=2,
                colour="red"
                ) +
                geom_text(data=data.frame(x=data@doseGrid,
                                          y=as.numeric(rr) 
                ), 
                mapping=aes(x=data@doseGrid, 
                            y=as.numeric(rr) , 
                            label=paste(round(as.numeric(rr), 1),"%")
                ), 
                size = 3,
                hjust = 0.5,
                vjust = -1
                )+
                geom_vline(xintercept = doseMTD,
                           linetype = "dashed",
                           size=1,
                           colour="blue"
                           ) +
                annotate(geom = "text", 
                         x = (doseMTD),
                         y = 97.5, 
                         label = "Max",
                         colour="blue"
                        ) +
                geom_vline(xintercept = ret,
                           linetype = "dashed",
                           size=1,
                           colour="green") +
                annotate(geom = "text", 
                         x = ret,
                         y = 92.5, 
                         colour="green",
                         label = "Next"
                ) +
                xlab("Dose") +
                ylab(paste("Allocation criterion [%]")) +
                scale_x_continuous(limits=range(data@doseGrid) , breaks=data@doseGrid) +
                scale_y_continuous(limits=c(0, 100) , breaks=seq(0,100,10))
              
              ## if doselimit exists add the dose range on the plot
              if(length(doselimit))
              {
                plot1 <- plot1+
                  geom_rect(aes(xmin = data@doseGrid[1], 
                                ymin = 0, 
                                xmax =doselimit, 
                                ymax = 100
                  ),
                  alpha = 0.5
                  ) + 
                  annotate(geom = "text", 
                           x = doselimit,
                           y = 95, 
                           label = "Doselimit"
                  ) 
              }
              
              ## return next best dose
              return(list(value=ret,
                          plot=plot1))
              
            })

