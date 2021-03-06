% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Model-class.R
\docType{class}
\name{ModelTox-class}
\alias{ModelTox-class}
\alias{.ModelTox}
\title{No intialization function
Class for DLE models using pseudo data prior.
This is a class of DLE (dose-limiting events) models/ toxicity model which contains all DLE models
for which their prior are specified in form of pseudo data (as if there is some data before
the trial starts). It inherits all slots from \code{\linkS4class{ModelPseudo}}}
\description{
The \code{data} must obey the convention that the data input is called exactly in the
\code{\linkS4class{Data}} class. This refers to any observed DLE responses (\code{y} in
\code{\linkS4class{Data}} class), the dose (levels) (\code{x} in \code{\linkS4class{Data}} class)
at which these responses are observed, all dose levels considered in the study (\code{doseGrid}
in \code{\linkS4class{Data}}) class and other specifications in \code{\linkS4class{Data}}
class that can be used to generate prior or
posterior modal estimates or samples estimates for model parameter(s). If no responses is observed,
at least \code{doseGrid} in \code{\linkS4class{Data}} has to be specified in \code{data} slot for which
prior modal estimates or samples can be obtained for model parameters based on the specified pseudo
data.
}
\section{Slots}{

\describe{
\item{\code{data}}{refers to the data input specification in \code{\linkS4class{Data}} class which are used to
obtain model parameters estimates or samples (see details above)}
}}

\seealso{
\code{\linkS4class{LogisticIndepBeta}},
\code{\linkS4class{Effloglog}},
\code{\linkS4class{EffFlexi}}
}
\keyword{classes}
