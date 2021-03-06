% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FGS_Normal.R
\name{FGS_Normal}
\alias{FGS_Normal}
\title{SSR for flexible group sequential multi-arm design with normal outcomes}
\usage{
FGS_Normal(a, sel_rule, n1 = NULL, N2 = NULL, initialSS = 2, r21, rmax,
  alpha = 0.025, ushape = "obf", lfix = 0, theta, sigma, x0_1, x1, beta,
  ConditionalPower, ssr, delta = NULL, side = "upper")
}
\arguments{
\item{a}{Number of experimental treatments (a = 2 or 3)}

\item{sel_rule}{Treatment selection rule at interim analysis. There are four treatment selection rule, 1) sel_tr = 1, select the best treatment; 2) sel_tr=2, keep all promising treatments that do not cross either efficacy boundaries or futility boundaries; 3) sel_tr = 3, select the best treatment only if the mean difference to control is above a threshold delta; 4) sel_tr = 4, keep all promising treatments where the mean difference to control is above delta.}

\item{n1}{First stage sample size per arm (default is NULL).}

\item{N2}{Total sample size per arm at the end of the trial (default is NULL).}

\item{initialSS}{Whether initial sample size needs to be calculated. If initial sample size does not need to be calculated, then n1 and N2 have to be provided and initialSS is set to 0. If initial sample size need to be calculated, then initialSS can take numeric value 1 (initial sample size will be calculated using the method for traditional fixed two-arm design) or 2 (initial sample size will be calculated using the method for group sequential two-arm design).}

\item{r21}{Allocation ratio of N2 to n1 (r21 = N2/n1). If initialSS = NULL, then it should also be NULL. If initialSS = 1 or 2, then it should be a numeric value that is larger than 1.}

\item{rmax}{The ratio of the maximum sample size allowed per arm (N2_max) to N2 (rmax = N2_max/N2). It should be a numeric value that is larger than 1.}

\item{alpha}{One sided familywise error rate (default is 0.025).}

\item{ushape}{The shape of the efficacy boundaries (default is "obf"). Either a function specifying the shape or one of "pocock", "obf" (the default), "triangular" and "fixed".}

\item{lfix}{Fixed lower boundary on z-scale (default is 0). lfix must be smaller than Phi(1-alpha)/2 to ensure that it is smaller than the upper boundary.}

\item{theta}{Planned absolute mean difference of experimental treatment to control. It is used for initial sample size caculation and SSR based on the planned treatment effect.}

\item{sigma}{Known standard deviation of normal distribution.}

\item{x0_1}{First stage sample mean of control arm. It has to be a one column matrix with the same number of rows as x1.}

\item{x1}{First stage sample mean of experimental treatment arms. It has to be a matrix with each row represents one trial, and the number of colum should be equal to a.}

\item{beta}{Type II error. 1-beta is used as the desired power and conditional power.}

\item{ConditionalPower}{Methods for conditional power. If ConditionalPower = "planned", then conditional power will be calculated using the planned treatment effect (theta). If ConditionalPower = "observed", then conditional power is calcuated using the maximum observed treatment effect at interim analysis.}

\item{ssr}{Whether sample size re-estimation should be performed. If ssr = 0, sample size re-estimation will not be performed. If ssr = 1, sample size will be performed if conditional power is lower than 1-beta.}

\item{delta}{The minimum treatment effect that is of interest. The default is NULL, but when sel_rule = 3 or 4, a numeric value has to be specified.}

\item{side}{Side for conditional power (default is "upper"). Either "upper" or "lower".}
}
\value{
A list containing the following components:
\item{EfficacyBound}{The efficacy boundaries for both interim analysis and final analysis.}
\item{FutilityBound}{The futility boundaries for both interim analysis and final analysis.}
\item{stop_efficacy}{A matrix indicating which arms are stopped for efficacy. 1 indicates stopping for efficacy and 0 otherwise.}
\item{stop_futility}{A matrix indicating which arms are stopped for futility. 1 indicates stopping for efficacy and 0 otherwise.}
\item{sel_tr}{A matrix indicating which experimental treatments can continue to the second stage.}
\item{cp}{Conditional power at interim analysis. If a trial stops at interim analysis, conditional power will not be calculated and returned as NA.}
\item{n1}{First stage sample size per arm.}
\item{n2new}{New second stage sample size per arm.}
\item{cp_new}{Conditional power based on the new second stage sample size.}
\item{NoRoot}{If new sample size is successful found, then NoRoot = NULL. If new sample size cann't be found, then NoRoot will return the index of the trial for which the new sample size is not found.}
}
\description{
Calculate conditional power and perform SSR at the interim analysis for flexible group sequential multi-arm design with normal outcomes. If conditional power is less than desired, perform SSR and calculate the new conditional power based on the new sample size.
}
\details{
This function calculates conditional power and perform sample size re-estimation (SSR) based on the conditional power approach for two-stage flexible group sequential multi-arm study design (as described in Magirr et al (2014)) with normal outcomes. It is assumed that sample size is equally allocated among arms and the observations from each arm follow a normal distribution with a common known variance. At interim analysis, a trial will be stopped if any treatments cross the efficacy boundaries and any treatments cross the futility boundary will be dropped. Conditional power is calculated under the global alternative hypothesis where it is assumed all the treatment effects are equal to planned treatment effect (theta) or maximum observed treatment effect. If a trial stops at interim analysis, conditional power(CP), new second stage sample size per arm (n2new) and new conditional power(cp_new) will not be calculated and will all be returned as NA. If a trial can continue to the second stage and conditional power is larger than 1-beta, then n2new and cp_new will not be calculated and will be returned as NA. Multiple trials can be assessed at the same time.
}
\examples{
two-stage three experiments versus control design, initial sample size was calculated based on group sequential design and the best treatment is selected at the interim analysis
FGS_Normal(a = 3, sel_rule = 2, n1=NULL, N2=NULL, initialSS = 2, r21 = 2, rmax= 1.5, alpha = 0.025, ushape = "obf", lfix=0, theta=2, sigma = 6, x0_1=matrix(2.2,nrow=2,ncol=1), x1=matrix(c(2,2.2,3.1,3.0,3.5,2.5), nrow=2,ncol=3), beta=0.2, ConditionalPower = "observed",ssr = 1,delta = NULL)
}
\references{
Magirr D, Jaki T, Whitehead J: A generalized Dunnett test for multi-arm multi-stage clinical studies with treatment selection. Biometrika 2012, 99(2):494-501.

Magirr D, Stallard N, Jaki T: Flexible sequential designs for multi-arm clinical trials. Statistics in Medicine 2014, 33(19):3269-3279.
}
