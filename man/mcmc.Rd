% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc.R
\name{mcmc}
\alias{mcmc}
\alias{mcmc,GeneralData,GeneralModel,McmcOptions-method}
\alias{mcmc-GeneralData}
\alias{mcmc,DataMixture,GeneralModel,McmcOptions-method}
\alias{mcmc,Data,LogisticIndepBeta,McmcOptions-method}
\alias{mcmc,DataDual,Effloglog,McmcOptions-method}
\alias{mcmc,DataDual,EffFlexi,McmcOptions-method}
\title{Obtaining Posterior Samples for all Model Parameters}
\usage{
mcmc(data, model, options, ...)

\S4method{mcmc}{GeneralData,GeneralModel,McmcOptions}(data, model, options, from_prior = data@nObs == 0L, ...)

\S4method{mcmc}{DataMixture,GeneralModel,McmcOptions}(
  data,
  model,
  options,
  from_prior = data@nObs == 0L & data@nObsshare == 0L,
  ...
)

\S4method{mcmc}{Data,LogisticIndepBeta,McmcOptions}(data, model, options, ...)

\S4method{mcmc}{DataDual,Effloglog,McmcOptions}(data, model, options, ...)

\S4method{mcmc}{DataDual,EffFlexi,McmcOptions}(data, model, options, ...)
}
\arguments{
\item{data}{(\code{GeneralData})\cr an input data.}

\item{model}{(\code{GeneralModel})\cr an input model.}

\item{options}{(\code{McmcOptions})\cr MCMC options.}

\item{...}{not used.}

\item{from_prior}{(\code{flag})\cr sample from the prior only? Default to \code{TRUE}
when number of observations in \code{data} is \code{0}. For some models it might be
necessary to specify it manually here though.}
}
\value{
The posterior samples, an object of class \code{\link{Samples}}.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

This is the function that actually runs the \code{JAGS} MCMC machinery to produce
posterior samples from all model parameters and required derived values.
It is a generic function, so that customized versions may be conveniently
defined for specific subclasses of \code{\link{GeneralData}}, \code{\link{GeneralModel}}, and
\code{\link{McmcOptions}} input.
}
\section{Functions}{
\itemize{
\item \code{mcmc,GeneralData,GeneralModel,McmcOptions-method}: Standard method which uses JAGS.

\item \code{mcmc,DataMixture,GeneralModel,McmcOptions-method}: Method for DataMixture with different from_prior default

\item \code{mcmc,Data,LogisticIndepBeta,McmcOptions-method}: Obtain posterior samples for the model parameters based on the pseudo 'LogisticsIndepBeta'
DLE model. The joint prior and posterior probability density function of
the intercept \eqn{\phi_1} (phi1) and the slope \eqn{\phi_2} (phi2) are given in Whitehead and
Williamson (1998) and TsuTakawa (1975). However, since asymptotically, the joint posterior probability density
will be bivariate normal and we will use the bivariate normal distribution to
generate posterior samples of the intercept and the slope parameters. For the prior samples of
of the intercept and the slope a bivariate normal distribution with mean and the covariance matrix given in Whitehead and
Williamson (1998) is used.

\item \code{mcmc,DataDual,Effloglog,McmcOptions-method}: Obtain the posterior samples for the model parameters in the
Efficacy log log model. Given the value of \eqn{\nu}, the precision of the efficacy responses,
the joint prior or the posterior probability of the intercept \eqn{\theta_1} (theta1) and
the slope \eqn{\theta_2} (theta2) is a bivariate normal distribution. The  \eqn{\nu} (nu),
the precision of the efficacy responses is either a fixed value or has a gamma distribution.
If a gamma distribution is used, the samples of nu will be first generated.
Then the mean of the of the nu samples
will be used the generate samples of the intercept and slope parameters of the model

\item \code{mcmc,DataDual,EffFlexi,McmcOptions-method}: Obtain the posterior samples for the estimates in the Efficacy Flexible form.
This is the mcmc procedure based on what is described in Lang and Brezger (2004) such that
samples of the mean efficacy responses at all dose levels, samples of sigma2 \eqn{sigma^2},
the variance of the efficacy response and samples of sigma2betaW \eqn{sigma^2_{beta_W}}, the variance of
the random walk model will
be generated. Please refer to Lang and Brezger (2004) for the procedures and the form of
the joint prior and posterior probability density for the mean efficacy responses. In addition,
both sigma2 and sigma2betaW can be fixed or having an inverse-gamma prior and posterior distribution.
Therefore, if the inverse gamma distribution(s) are used, the parameters in the distribution will be
first updated and then samples of sigma2 and sigma2betaW will be generated using the updated parameters.
}}

\note{
The type of Random Number Generator (RNG) and its initial seed used by
\code{JAGS} are taken from the \code{options} argument. If no initial values are
supplied (i.e RNG kind or seed slot in \code{options} has \code{NA}), then they will
be generated automatically by \code{JAGS}.
}
\examples{
# Create some data from the class `Data`.
my_data <- Data(
  x = c(0.1, 0.5, 1.5, 3, 6, 10, 10, 10),
  y = c(0, 0, 0, 0, 0, 0, 1, 0),
  doseGrid = c(0.1, 0.5, 1.5, 3, 6, seq(from = 10, to = 80, by = 2))
)

# Initialize the CRM model.
my_model <- LogisticLogNormal(
  mean = c(-0.85, 1),
  cov = matrix(c(1, -0.5, -0.5, 1), nrow = 2),
  ref_dose = 56
)

# Sample from the posterior distribution.
my_options <- McmcOptions(
  burnin = 100,
  step = 2,
  samples = 1000
)

samples <- mcmc(data = my_data, model = my_model, options = my_options)
samples
##obtain mcmc DLE samples given the data, LogisticIndepBeta (DLE model) and mcmc simulations options
## data must be of 'Data' class
data<-Data(x=c(25,50,50,75,100,100,225,300),y=c(0,0,0,0,1,1,1,1),
           doseGrid=seq(25,300,25))
## model must be of 'LogisticIndepBeta' class
model<-LogisticIndepBeta(binDLE=c(1.05,1.8),DLEweights=c(3,3),DLEdose=c(25,300),data=data)
## options must be ''McmcOptions' class
options<-McmcOptions(burnin=100,step=2,samples=200)
set.seed(94)
samples<-mcmc(data=data,model=model,options=options)
##obtain mcmc efficacy samples given the data, 'Effloglog' model (efficacy model) and
## mcmc simulations options data must be of 'DataDual' class
data<-DataDual(x=c(25,50,25,50,75,300,250,150),
              y=c(0,0,0,0,0,1,1,0),
              w=c(0.31,0.42,0.59,0.45,0.6,0.7,0.6,0.52),
              doseGrid=seq(25,300,25),placebo=FALSE)
## model must be of 'Effloglog' class
Effmodel<-Effloglog(Eff=c(1.223,2.513),Effdose=c(25,300),nu=c(a=1,b=0.025),data=data,c=0)

## options must be ''McmcOptions' class
options<-McmcOptions(burnin=100,step=2,samples=200)
set.seed(94)
samples<-mcmc(data=data,model=Effmodel,options=options)
##obtain mcmc efficacy samples given the data, 'EffFlexi' model (efficacy model) and 
## mcmc simulations options
## data must be of 'DataDual' class
data<-DataDual(x=c(25,50,25,50,75,300,250,150),
               y=c(0,0,0,0,0,1,1,0),
               w=c(0.31,0.42,0.59,0.45,0.6,0.7,0.6,0.52),
               doseGrid=seq(25,300,25))
## model must be of 'EffFlexi' class

Effmodel<- EffFlexi(Eff=c(1.223, 2.513),Effdose=c(25,300),
                    sigma2=c(a=0.1,b=0.1),sigma2betaW=c(a=20,b=50),smooth="RW2",data=data)

## options must be ''McmcOptions' class
options<-McmcOptions(burnin=100,step=2,samples=200)
set.seed(94)
samples<-mcmc(data=data,model=Effmodel,options=options)
}
