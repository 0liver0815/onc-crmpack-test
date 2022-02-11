# You can install the development version of crmPack from github with:
# install.packages("devtools")
# devtools::install_github("Roche/crmPack")

#install and subsequently load package in R
# You can install the stable release version of crmPack from CRAN with:
# install.package("crmPack")
library(crmPack)

# install.packages("ggmcmc") to be used for plotting results produced by crmPack
library(ggmcmc)

###### NOTE ####################################################################################################
## open and run the file helpers.R in the in the crmPack installation R folder, before proceeding further below
## this is file contains Some helper functions used during development of new S4 class objects further below
###### NOTE ####################################################################################################

# setwd("/home/crm/dk/dev")
# getwd() 
# 
# source("LogisticKadane2.R")





emptydata <- Data(doseGrid = c(1.5, 2.5, 3.5, 4.5, 6, 7 ))

MSLN_TTC_Model <- LogisticKadane2(theta = 0.3, 
                                      dmin = 1.5, 
                                      dmax = 7, 
                                      alpha = 1, 
                                      beta = 19, 
                                      shape = 0.5625, 
                                      rate = 0.125)

options <- McmcOptions(burnin = 1000, step = 1, samples = 1000)
set.seed(94)

# Prior exploration 
PriorSamples <- mcmc(emptydata, MSLN_TTC_Model, options)

plot(PriorSamples, MSLN_TTC_Model, emptydata) + ggtitle("Prior")

p0samples <- get(PriorSamples, "p0")
MTDsamples <- get(PriorSamples, "MTD")


print(ggs_traceplot(p0samples)+ ggtitle("Prior"))
print(ggs_traceplot(MTDsamples)+ ggtitle("Prior"))

print(ggs_autocorrelation(p0samples)+ ggtitle("Prior"))
print(ggs_autocorrelation(MTDsamples)+ ggtitle("Prior"))