simsum <- function(scenario, object, truth, target, trueMTD){
## extract dose grid
doseGrid <- object@data[[1]]@doseGrid

## evaluate true toxicity at doseGrid
trueTox <- truth(doseGrid)

## what are the levels above target interval?
xAboveTarget <- which(trueTox > target[1])

##number of simulations
nsims <- length(object@doses)

## MTD and CV
MTD <- sapply(object@fit, function(x) {x[length(doseGrid)+1,"middle"]})
CV <- sapply(object@fit, function(x) {x[length(doseGrid)+1,"lower"]})

## Stop reasons
#split stop reason in two parts
stop <- strsplit(unlist(object@stopReasons),' : ')
#print(stop[[1:5]])
#convert into logical vector
stop1 <- as.logical(lapply(stop, function(l) l[[1]]))
#put vector into a matrix
stop.res <- matrix(stop1,nrow=length(object@doses),byrow=T)
#calculate the stop reasons for each column of the matrix
#several stop reasons may occur at the same time
#the columns appear in the order as they are out together in the simulation program
#myStoppinglow | myStoppinghigh | myStoppingCV | myStoppingnpat | myStoppingfirst
stopreason <- apply(stop.res,2,sum)/nsims*100
names(stopreason) <- c('tox rule', 'tox mod', 'saf mod', 'CV', 'n9')

## count only the first occurence, i.e. stop reason needs to be in the correct order
stop.res2 <- apply(stop.res,1,which)
stop2 <- sapply(stop.res2, function(l) l[[1]])
stop2 <- factor(stop2,levels=1:5)
stopreason2 <- table(stop2)/nsims*100
row.names(stopreason2) <- c('tox rule', 'tox mod', 'saf mod', 'CV', 'n9')

##Trials where MTD estimate is below first or above top dose
above.dmax <- mean(MTD>tail(doseGrid, n=1))*100
below.dmin <- mean(MTD<doseGrid[1])*100

## MTD: median, rel error, CV
#mean(MTD)
MTD.med <- median(MTD)
#median relative error
MTD.relerr <- median(abs(trueMTD-MTD)/trueMTD)*100
MTD.CV <- mean(CV)*100

## number of patients overdosed
nAboveTarget <- sapply(object@data,
                       function(d){
                         sum(d@xLevel %in% xAboveTarget)
                       })
#nAboveTarget
od.n <- mean(nAboveTarget)
od.sd <- sd(nAboveTarget)
#min(nAboveTarget)
#max(nAboveTarget)
od.range <- range(nAboveTarget)

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

return(c(scenario=scenario,
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
         n=n,
         n.sd=n.sd,
         n.range=n.range))
}

load("~/results/18795_msln_results_seed9.RData")


#simsum(1, safe.18795, safe, 0.3, 15)

res_seed9 <-
t(
rbind(
  simsum(1, safe.18795, safe, 0.3, 15),
  simsum(2, late.18795, late, 0.3, 4.5),
  simsum(3, early.18795, early, 0.3, 2.5),
  simsum(4, toxic.18795, toxic, 0.3, 1),
  simsum(5, peak.18795, peak, 0.3, 3.5)
))

#format(round(res,2),nsmal=2)
format(round(res_seed9,2),nsmal=2)

# plot the course of the 5th simulated trial as follows:
for (i in 1:5){print(plot(safe.18795@data[[i]]))}
for (i in 1:5){print(plot(late.18795@data[[i]]))}
for (i in 1:5){print(plot(early.18795@data[[i]]))}
for (i in 1:5){print(plot(toxic.18795@data[[i]]))}
for (i in 1:5){print(plot(peak.18795@data[[i]]))}

# The final dose for this trial was
safe.18795@doses[1:5]

# and the stopping reason was
safe.18795@stopReasons[[5]]

# plot all simulations for this scenario
# The resulting plot shows on the top panel a summary of the trial trajectories. On the
# bottom, the proportions of doses tried, averaged over the simulated trials, are shown
print(plot(safe.18795))

simSum <- summary(safe.18795, truth=safe)

simSum

print(plot(simSum))

dosePlot <- plot(simSum, type="doseSelected") +
  scale_x_continuous(breaks=1:8, limits=c(0, 10))

print(dosePlot)

mean( unlist(safe.18795@fit[[1]]) )



    
