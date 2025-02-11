% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GDILM_SEIRS_Sim_Par_Est.R
\name{GDILM_SEIRS_Sim_Par_Est}
\alias{GDILM_SEIRS_Sim_Par_Est}
\title{GDILM SEIRS for a Simulation Study}
\usage{
GDILM_SEIRS_Sim_Par_Est(
  GridDim1,
  GridDim2,
  NPostPerGrid,
  MaxTimePand,
  tau0,
  lambda0,
  alphaS0,
  delta0,
  alphaT0,
  PopMin,
  PopMax,
  InfFraction,
  ReInfFraction,
  InfPrd,
  IncPrd,
  NIterMC,
  NIterMCECM
)
}
\arguments{
\item{GridDim1}{First dimension of the grid}

\item{GridDim2}{Second dimension of the grid}

\item{NPostPerGrid}{Number of postal codes per grid cell}

\item{MaxTimePand}{Last time point of the pandemic}

\item{tau0}{Initial value for spatial precision}

\item{lambda0}{Initial value for spatial dependence}

\item{alphaS0}{Initial value for the susceptibility intercept}

\item{delta0}{Initial value for the spatial decay parameter}

\item{alphaT0}{Initial value for the infectivity intercept}

\item{PopMin}{Minimum population per postal code}

\item{PopMax}{Maximum population per postal code}

\item{InfFraction}{Fraction of each grid cell's population to be infected}

\item{ReInfFraction}{Fraction of each grid cell's population to be reinfected}

\item{InfPrd}{Infectious period that can be obtained either from the literature or by fitting an SEIRS model to the data}

\item{IncPrd}{Incubation period that can be obtained either from the literature or by fitting an SEIRS model to the data}

\item{NIterMC}{Number of MCMC iterations}

\item{NIterMCECM}{Number of MCECM iterations}
}
\value{
\code{alphaS} Estimate of alpha S

\code{BetaCovInf} Estimate of beta vector for the individual level infection covariate

\code{BetaCovSus} Estimate of beta vector for the areal susceptibility to first infection covariate

\code{BetaCovSusReInf} Estimate of beta vector for the areal susceptibility to reinfection covariate

\code{alphaT} Estimate of alpha T

\code{delta} Estimate of delta

\code{tau1} Estimate of tau

\code{lambda1} Estimate of lambda

\code{AIC} AIC of the fitted GDILM SEIRS
}
\description{
This function conducts a simulation study for the Geographically Dependent Individual Level Model (GDILM) of infectious disease transmission, incorporating reinfection dynamics within the Susceptible-Exposed-Infectious-Recovered-Susceptible (SEIRS) framework, using a user-defined grid size. It applies a likelihood based Monte Carlo Expectation Conditional Maximization (MCECM) algorithm to estimate model parameters and compute the AIC.
}
\examples{
\donttest{
GDILM_SEIRS_Sim_Par_Est(5,5,10,30,0.7, 0.5, -1, 2.5, 0,40, 50,0.3,0.6, 3, 3, 10, 3)
}

}
