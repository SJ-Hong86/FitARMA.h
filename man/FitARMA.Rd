% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FitARMA.R
\name{FitARMA}
\alias{FitARMA}
\title{Fit ARMA/ARIMA using fast MLE algorithm}
\usage{
FitARMA(z, order = c(0, 0, 0), demean = TRUE, 
  MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
}
\arguments{
\item{z}{time series.}

\item{order}{model order, c(p,d,q).}

\item{demean}{if TRUE, mean parameter included otherwise assumed zero.}

\item{MeanMLEQ}{exact MLE for mean, ignored unless demean=TRUE.}

\item{pApprox}{order of approximation to be used.}

\item{MaxLag}{maximum number of lags for portmanteau test.}
}
\value{
A list with class name "FitARMA" and components:
\itemize{
\item \code{loglikelihood} value of the loglikelihood.
\item \code{phiHat} AR coefficients.
\item \code{thetaHat} MA coefficients.
\item \code{sigsqHat} innovation variance estimate.
\item \code{muHat} estimate of the mean.
\item \code{covHat} covariance matrix of the coefficient estimates.
\item \code{racf} residual autocorrelations.
\item \code{LjungBox} table of Ljung-Box portmanteau test statistics.
\item \code{res} innovation residuals, same length as z.
\item \code{fits} fitted values, same length as z.
\item \code{demean} TRUE if mean estimated otherwise assumed zero.
\item \code{IterationCount} number of iterations in mean mle estimation.
\item \code{convergence} value returned by optim – should be 0.
\item \code{MLEMeanQ} TRUE if mle for mean algorithm used.
\item \code{tsp} tsp(z).
\item \code{call} result from match.call() showing how the function was called.
\item \code{ModelTitle} description of model.
\item \code{DataTitle} returns attr(z,"title").
}
}
\description{
Fits an ARIMA(p,d,q) model using the algorithm given in McLeod and Zhang (2007).
}
\details{
See McLeod and Ying (2007).
}
\note{
When d>0 and demean=TRUE, the mean of the differenced series is estimated.
This corresponds to including a polynomial of degree d.

When d>0, the AIC/BIC are computed for the differenced series and so
they are not comparable to the values obtained for models with d=0.
}
\examples{
data(SeriesA) #in datasets()
out1<-FitARMA(SeriesA, c(1,0,1))
out1
coef(out1)
out2<-FitARMA(SeriesA, c(0,1,1))
out2
coef(out2) 

}
\references{
A.I. McLeod andY. Zhang (2008), Faster ARMA maximum likelihood estimation,
Computational Statistics & Data Analysis, 52-4, 2166-2176. DOI link:
http://dx.doi.org/10.1016/j.csda.2007.07.020.
}
\seealso{
\code{\link[=GetFitARMA]{GetFitARMA()}}, \code{\link[=print.FitARMA]{print.FitARMA()}}, \code{\link[=coef.FitARMA]{coef.FitARMA()}}, \code{\link[=residuals.FitARMA]{residuals.FitARMA()}},
\code{\link[=fitted.FitARMA]{fitted.FitARMA()}}, \code{\link[=arima]{arima()}}.
}
\author{
A.I. McLeod.
}
