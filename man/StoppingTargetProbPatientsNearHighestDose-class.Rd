% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/ZZ_StoppingTargetProbPatientsNearHighestDose_DK.R
\docType{class}
\name{StoppingTargetProbPatientsNearHighestDose-class}
\alias{StoppingTargetProbPatientsNearHighestDose-class}
\alias{.StoppingTargetProbPatientsNearHighestDose}
\title{Stop based on number of patients near to next best dose}
\description{
Stop based on number of patients near to next best dose
}
\section{Slots}{

\describe{
\item{\code{target}}{the target toxicity interval, e.g. \code{c(0.2, 0.35)}}

\item{\code{prob}}{required target toxicity probability (e.g. \code{0.4})
for reaching sufficient precision}

\item{\code{nPatients}}{number of required patients}

\item{\code{percentage}}{percentage (between 0 and 100) within the next best dose
the patients must lie}
}}

\keyword{classes}
