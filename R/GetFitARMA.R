#' Fit ARMA(p,q) model to mean zero time series.
#' 
#' The algorithm of McLeod and Zhang (2007) is used.
#' 
#' @usage GetFitARMA(y, p, q, pApprox = 30, init = 0)
#' @param y time series.
#' @param p AR order.
#' @param q MA order.
#' @param pApprox AR approximation.
#' @param init initial parameter estimates.
#' @details See McLeod and Zhang (2006).
#' @returns
#' A list with class name "FitARMA" and components:
#' 
#' * `loglikelihood` value of maximized loglikelihood.
#' * `phiHat` estimated phi parameters.
#' * `thetaHat` estimated theta parameters.
#' * `convergence` result from optim.
#' * `algorithm` indicates "L-BFGS-B" or "Nelder-Mead" according as which algorithm was used in optim.
#' @author A.I. McLeod.
#' @references
#' A.I. McLeod andY. Zhang (2008), Faster ARMA maximum likelihood estimation, 
#' Computational Statistics & Data Analysis, 52-4, 2166-2176. DOI link: 
#' http://dx.doi.org/10.1016/j.csda.2007.07.020.
#' @seealso [arima()], [FitARMA()].
#' @examples 
#' data(SeriesA)
#' z<-SeriesA-mean(SeriesA)
#' GetFitARMA(z, 1, 1)
#' w<-diff(z, differences=1)
#' GetFitARMA(w, 0, 1) 
#' 
#' @export
GetFitARMA <-
  function(y, p, q, pApprox=30, init=0){
    #White noise case
    if (p==0 && q==0) {
      LogL<-LoglikelihoodAR(0, y)
      ans<-list(loglikelihood=LogL, phiHat=numeric(0), thetaHat=numeric(0),convergence=NULL,algorithm=NULL)
      return(ans)
    }
    if (length(init)==1 && init == 0) 
      xinit<-numeric(p+q) #handles both p+q>1 and p+q=1
    else
      xinit<-init
    if (length(xinit)!=p+q){
      message("GetARMAFit: init length not correct. Set to zero.")
      xinit<-numeric(p+q)
    }
    if (max(abs(xinit))>0.99){
      message("GetARMAFit: init parameter setting outside (-0.99,0.99). Reset to 0.")
      xinit<-numeric(p+q)
    }
    alg<-1 
    n<-length(y)  
    penaltyLoglikelihood<-(-length(y)/2*log(sum(y^2)/n))-10^4 
    #setup Champernowne matrix
    if (q>0)
      pApp <- pApprox
    else
      pApp <- p
    CD<-ChampernowneD(y,pApp,MeanZero=TRUE)
    n<-length(y)
    
    #ARMA case
    if (p>0 && q>0) {
      EntropyARMA<-function(x){
        if (max(abs(x))>0.99) 
          -penaltyLoglikelihood+max(abs(x))^2 
        else {
          zetaPhi=x[1:p]
          zetaTheta=x[(p+1):(p+q)]
          g<-TacvfARMA(PacfToAR(zetaPhi),PacfToAR(zetaTheta),pApp)
          xpar<-PacfDL(g, LinearPredictor=TRUE)$ARCoefficients
          -FastLoglikelihoodAR(xpar,n,CD)
        }
      }
      ans<-optim(xinit,EntropyARMA,method="L-BFGS-B", lower=rep(-0.99,p+q), upper=rep(0.99,p+q),
                 control=list(trace=0))
      if(ans$convergence != 0  || max(abs(ans$par))>0.99) {#convergence problem. Use Nelder-Mead with penalty function
        alg<-2
        ans<-optim(c(0,0),EntropyARMA,method="Nelder-Mead")
      }
      zpar<-ans$par
      phiHat<-PacfToAR(zpar[1:p])
      thetaHat<-PacfToAR(zpar[p+(1:q)])
    }
    
    #AR Case
    if (p>0 && q==0) {
      EntropyARMA<-function(x){
        if (max(abs(x))>0.99) 
          -penaltyLoglikelihood 
        else {
          zetaPhi=x[1:p]
          xpar<-PacfToAR(zetaPhi)
          -FastLoglikelihoodAR(xpar,n,CD)
        }
      }
      ans<-optim(xinit,EntropyARMA,method="L-BFGS-B", lower=rep(-0.99,p+q), upper=rep(0.99,p+q))
      if(ans$convergence != 0) {
        alg<-2
        ans<-optim(xinit,EntropyARMA,method="Nelder-Mead")
      }
      zpar<-ans$par
      phiHat<-PacfToAR(zpar[1:p])
      thetaHat<-numeric(0)
    }
    
    #MA Case
    if (p==0 && q>0) {
      EntropyARMA<-function(x){
        if (max(abs(x))>0.99) 
          -penaltyLoglikelihood 
        else {
          zetaTheta=x[1:q]
          g<-TacvfARMA(numeric(0),PacfToAR(zetaTheta),pApp)
          xpar<-PacfDL(g, LinearPredictor=TRUE)$ARCoefficients
          -FastLoglikelihoodAR(xpar,n,CD)
        }
      }
      ans<-optim(xinit,EntropyARMA,method="L-BFGS-B", lower=rep(-0.99,p+q), upper=rep(0.99,p+q))
      if(ans$convergence != 0) {
        alg<-2
        ans<-optim(xinit,EntropyARMA,method="Nelder-Mead")
      }
      zpar<-ans$par
      phiHat<-numeric(0)
      thetaHat<-PacfToAR(zpar[1:q])
    }
    list(loglikelihood=-ans$value, phiHat=phiHat, thetaHat=thetaHat, convergence=ans$convergence, 
         algorithm=c("L-BFGS-B","Nelder-Mead")[alg])
  }
