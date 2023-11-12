#' Fit ARMA/ARIMA using fast MLE algorithm
#' 
#' Fits an ARIMA(p,d,q) model using the algorithm given in McLeod and Zhang (2007).
#' 
#' @usage FitARMA(z, order = c(0, 0, 0), demean = TRUE, 
#'   MeanMLEQ = FALSE, pApprox = 30, MaxLag = 30)
#' @param z time series.
#' @param order model order, c(p,d,q).
#' @param demean if TRUE, mean parameter included otherwise assumed zero.
#' @param MeanMLEQ exact MLE for mean, ignored unless demean=TRUE.
#' @param pApprox order of approximation to be used.
#' @param MaxLag maximum number of lags for portmanteau test.
#' @details See McLeod and Ying (2007).
#' @returns
#' A list with class name "FitARMA" and components:
#' 
#' * `loglikelihood` value of the loglikelihood.
#' * `phiHat` AR coefficients.
#' * `thetaHat` MA coefficients.
#' * `sigsqHat` innovation variance estimate.
#' * `muHat` estimate of the mean.
#' * `covHat` covariance matrix of the coefficient estimates.
#' * `racf` residual autocorrelations.
#' * `LjungBox` table of Ljung-Box portmanteau test statistics.
#' * `res` innovation residuals, same length as z.
#' * `fits` fitted values, same length as z.
#' * `demean` TRUE if mean estimated otherwise assumed zero.
#' * `IterationCount` number of iterations in mean mle estimation.
#' * `convergence` value returned by optim – should be 0.
#' * `MLEMeanQ` TRUE if mle for mean algorithm used.
#' * `tsp` tsp(z).
#' * `call` result from match.call() showing how the function was called.
#' * `ModelTitle` description of model.
#' * `DataTitle` returns attr(z,"title").
#' @note 
#' When d>0 and demean=TRUE, the mean of the differenced series is estimated. 
#' This corresponds to including a polynomial of degree d.
#' 
#' When d>0, the AIC/BIC are computed for the differenced series and so 
#' they are not comparable to the values obtained for models with d=0.
#' @author A.I. McLeod.
#' @references
#' A.I. McLeod andY. Zhang (2008), Faster ARMA maximum likelihood estimation, 
#' Computational Statistics & Data Analysis, 52-4, 2166-2176. DOI link: 
#' http://dx.doi.org/10.1016/j.csda.2007.07.020.
#' @seealso [GetFitARMA()], [print.FitARMA()], [coef.FitARMA()], [residuals.FitARMA()],
#'   [fitted.FitARMA()], [arima()].
#' @examples 
#' data(SeriesA) #in datasets()
#' out1<-FitARMA(SeriesA, c(1,0,1))
#' out1
#' coef(out1)
#' out2<-FitARMA(SeriesA, c(0,1,1))
#' out2
#' coef(out2) 
#' 
#' @export
FitARMA <-
  function(z,order=c(0,0,0),demean=TRUE,MeanMLEQ=FALSE,pApprox=30,MaxLag=30){
    p<-order[1]
    d<-order[2]
    q<-order[3]
    Z<-z
    if (d > 0) Z<-diff(z, differences=d)
    if (demean) 
      mz<-mean(Z)
    else
      mz<-0
    y<-Z-mz
    pApp<-pApprox
    if (q == 0) pApp<-p
    ans<-GetFitARMA(y,p,q,pApp)
    LL<-ans$loglikelihood
    mu<-iter<-0
    if (MeanMLEQ && (p>0||q>0)) {
      etol <- MaxIter <- 10
      while(etol> 1e-06 && iter<MaxIter){
        LLPrev<-LL
        iter<-iter+1
        g<-TacvfARMA(ans$phi,ans$theta,pApp)
        coefAR<-PacfDL(g, LinearPredictor=TRUE)$ARCoefficients
        mu<-GetARMeanMLE(y, coefAR)        
        ans<-GetFitARMA(y-mu,p,q,pApp)
        LL<-ans$loglikelihood
        etol<-abs(LL-LLPrev)/LLPrev
        if (ans$convergence != 0) 
          stop("GetARFit returned convergence = ",ans$convergence)           
      }
    }
    muHat<-mu+mz
    phiHat<-ans$phiHat
    thetaHat<-ans$thetaHat
    if (p>0 || q>0) {
      g<-TacvfARMA(phiHat,thetaHat,pApp)
      coefAR<-PacfDL(g, LinearPredictor=TRUE)$ARCoefficients
      res<-BackcastResidualsAR(y, coefAR, Q=100, demean=FALSE)
    }
    else
      res<-Z-muHat
    fits<-Z-res
    n<-length(res)
    sigsq<-sum(res^2)/n
    if (p>0 || q>0) 
      covHat<-solve(InformationMatrixARMA(phiHat,thetaHat))/n
    else 
      covHat<-numeric(0)
    racf<-(acf(res, plot=FALSE, lag.max=MaxLag)$acf)[-1]
    LBQ<-LjungBoxTest(res, lag.max=MaxLag, k=order[1]+order[3])
    if (d == 0)
      if (q == 0) ModelTitle<-paste("AR(",p,")",sep="")
    if (q > 1) ModelTitle<-paste("ARMA(",p,",",q,")",sep="")
    else
      ModelTitle<-paste("ARIMA(",p,",",d,",",q,")",sep="")
    out<-list(loglikelihood=ans$loglikelihood,phiHat=phiHat,thetaHat=thetaHat,sigsqHat=sigsq,muHat=muHat,covHat=covHat,
              racf=racf, LjungBoxQ=LBQ,res=res,fits=fits, demean=demean,
              iterationCount=iter,convergence=ans$convergence, MeanMLE=MeanMLEQ, tsp=tsp(z),order=order,call=match.call(),
              DataTitle=attr(z,"title"),ModelTitle=ModelTitle)
    class(out)<-"FitARMA"
    out
  }
