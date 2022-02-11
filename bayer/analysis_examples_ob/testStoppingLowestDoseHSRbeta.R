#library(crmPack)

emptydata <- Data(doseGrid = c(1.5, 2.5, 3.5, 4.5, 6, 7))
emptydata.pl <- Data(doseGrid = c(0.001, 1.5, 2.5, 3.5, 4.5, 6, 7), placebo=T)

model <- LogisticKadaneBayer(theta = 0.3, 
                                      dmin = 1.5, 
                                      dmax = 7, 
                                      alpha = 1, 
                                      beta = 19, 
                                      shape = 0.5625, 
                                      rate = 0.125)

options <- McmcOptions(burnin = 10000, step = 2, samples = 10000)
set.seed(94)

# Prior exploration 
PriorSamples <- mcmc(emptydata, model, options)

plot(PriorSamples, MSLN_TTC_Model, emptydata) + ggtitle("Prior")

data0 <- update(emptydata, x=2.5, y=c(0,0,0))
data1 <- update(emptydata, x=1.5, y=c(0,0,1))
data2 <- update(emptydata, x=1.5, y=c(0,1,1))
data3 <- update(emptydata, x=1.5, y=c(1,1,1))

data.pl0 <- update(emptydata.pl, x=2.5, y=c(0,0,0))
data.pl2 <- update(emptydata, x=1.5, y=c(0,1,1))

stopping_o <- StoppingHardFirst(target=0.3, prob=0.9, label='test')
stopping_d <- StoppingLowestDoseHSRBeta_DK(target=0.3, prob=0.9, a=1, b=1)
stopping_n <- StoppingLowestDoseHSRBeta(target=0.3, prob=0.9) 

stopTrial(stopping_o, 7, PriorSamples, model, data0)
stopTrial(stopping_d, 7, PriorSamples, model, data0)
stopTrial(stopping_n, 7, PriorSamples, model, data0)

stopTrial(stopping_o, 7, PriorSamples, model, data1)
stopTrial(stopping_d, 7, PriorSamples, model, data1)
stopTrial(stopping_n, 7, PriorSamples, model, data1)

stopTrial(stopping_o, 7, PriorSamples, model, data2)
stopTrial(stopping_d, 7, PriorSamples, model, data2)
stopTrial(stopping_n, 7, PriorSamples, model, data2)

stopTrial(stopping_o, 7, PriorSamples, model, data3)
stopTrial(stopping_d, 7, PriorSamples, model, data3)
stopTrial(stopping_n, 7, PriorSamples, model, data3)

stopTrial(stopping_o, 7, PriorSamples, model, data.pl0)
stopTrial(stopping_d, 7, PriorSamples, model, data.pl0)
stopTrial(stopping_n, 7, PriorSamples, model, data.pl0)

stopTrial(stopping_o, 7, PriorSamples, model, data.pl2)
stopTrial(stopping_d, 7, PriorSamples, model, data.pl2)
stopTrial(stopping_n, 7, PriorSamples, model, data.pl2)
