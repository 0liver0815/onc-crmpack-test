#function that analyze the simulation results from crmpack

## for test purposes
# res.18795 <- readRDS(file = '~/test.RData')
# object=peak.18795
# truth=peak
# target=0.3
# trueMTD=3.5
# scenario=5

simsum <- function(object, scenario, truth, target, trueMTD){
  ## extract dose grid
  doseGrid <- object@data[[1]]@doseGrid

  ## evaluate true toxicity at doseGrid
  trueTox <- truth(doseGrid)

  ## closest to target
  closest <- min(which(abs(trueTox-target)==min(abs(trueTox-target))))
  if (trueTox[1] > target){closest = 0
  }else if (tail(trueTox,1) < target){closest = length(doseGrid)+1}

  ## what are the levels above target interval?
  ## important to round!!!!
  xAboveTarget <- which(round(trueTox,2) > target[1])

  ## number of simulations
  nsims <- length(object@doses)

  ## MTD and CV
  MTD <- sapply(object@estimates, function(l) l[[1]])
  CV <-  sapply(object@estimates, function(l) l[[2]])

  ## Stop reasons
  ## split stop reason in three parts
  stop <- strsplit(unlist(object@stopReasons),' : ')
  stop <- lapply(stop, trimws)
  #print(stop[[1]])
  ## convert into logical vector
  stop1 <- as.logical(lapply(stop, function(l) l[[1]]))

  ## extract label of stop reason
  stop.label <- strsplit(unlist(object@stopReasons[[1]]),' : ')

  ## put vector into a matrix
  stop.res <- matrix(stop1,nrow=length(object@doses),byrow=T)
  ## calculate the stop reasons for each column of the matrix
  ## several stop reasons may occur at the same time
  ## the columns appear in the order as they are out together in the simulation program
  stopreason <- apply(stop.res,2,sum)/nsims*100
  names(stopreason) <- sapply(stop.label, function(l) l[[2]])

  ## count only the first occurence, i.e. stop reason needs to be in the correct order
  stop.res2 <- apply(stop.res,1,which)
  stop2 <- sapply(stop.res2, function(l) l[[1]])
  stop2 <- factor(stop2,levels=1:length(stopreason))
  stopreason2 <- table(stop2)/nsims*100
  row.names(stopreason2) <- sapply(stop.label, function(l) l[[2]])

  ## Trials where MTD estimate is below first or above top dose
  above.dmax <- mean(MTD>tail(doseGrid, n=1))*100
  below.dmin <- mean(MTD<doseGrid[1])*100

  ## MTD: median, rel error, CV
  ## mean(MTD)
  MTD.med <- median(MTD)
  ## median relative error
  MTD.relerr <- median(abs(trueMTD-MTD)/trueMTD)*100
  MTD.CV <- mean(CV)*100

  ## number of patients overdosed
  nAboveTarget <- sapply(object@data,
                         function(d){
                           sum(d@xLevel %in% xAboveTarget)
                         })
  ## nAboveTarget
  od.n <- mean(nAboveTarget)
  od.sd <- sd(nAboveTarget)
  #min(nAboveTarget)
  #max(nAboveTarget)
  od.range <- range(nAboveTarget)
  od.quant <- quantile(nAboveTarget,c(0.01,0.99))
  #table(nAboveTarget)

  ## number of patients overall
  nObs <- sapply(object@data,
                 slot,
                 "nObs")
  #object@data[[1]]@nObs
  #nObs
  n <- mean(nObs)
  n.sd <- sd(nObs)
  #min(nObs)
  #max(nObs)
  n.range <- range(nObs)
  n.quant <- quantile(nObs,c(0.01,0.99))
  #table(nObs)

  ret <- c(scenario=scenario,
           stop1=stopreason,
           stop2=stopreason2,
           above.max=above.dmax,
           below.min=below.dmin,
           MTD.median=MTD.med,
           MTD.relerr=MTD.relerr,
           MTD.CV=MTD.CV,
           overd.n=od.n,
           overd.sd=od.sd,
           overd.range=od.range,
           overd.quant=od.quant,
           n=n,
           n.sd=n.sd,
           n.range=n.range,
           n.quant=n.quant)


  return(ret)
}
