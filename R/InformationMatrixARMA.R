#' Expected large-sample information matrix for ARMA
#' 
#' The expected large-sample information matrix per observation for ARMA(p,q) models is computed.
#' 
#' @usage InformationMatrixARMA(phi = numeric(0), theta = numeric(0))
#' @param phi AR coefficients.
#' @param theta MA coefficients.
#' @details The information matrix is derived by Box and Jenkins (1970).
#' @returns a matrix of order (p+q).
#' @author A.I. McLeod.
#' @references Box and Jenkins (1970). Time Series Analysis: Forecasting and Control.
#' @seealso [FitARMA()].
#' @examples 
#' # The covariance matrix estimates of the parameters phi and theta in an ARMA(1,1)
#' # with phi=0.9 and theta=0.5 and n=200 is
#' v<-solve(InformationMatrixARMA(0.9,0.5))/200
#' v
#' # and the standard errors are
#' sqrt(diag(v))
#' 
#' @export
InformationMatrixARMA <-
  function(phi = numeric(0), theta = numeric(0))
  {
    ########information matrix of arma#####################
    if(!(InvertibleQ(phi) & InvertibleQ(theta))) {
      cat("Model is non-causal or non-invertible\n")
      return(NULL)
    }
    p <- length(phi)
    q <- length(theta)
    if (p==0 && q==0) return(numeric(0))
    unames <- character(0)
    vnames <- character(0)
    if(p > 0) {
      if(p > 1) {
        vvmatrix <- (tccfAR(phi, phi)[ - (1:(p - 1))])[ - (p + 
                                                             1)]
      }
      else if(p == 1) {
        vvmatrix <- tccfAR(phi, phi)[ - (p + 1)]
      }
      vvmatrix <- toeplitz(vvmatrix)
      imatrix <- vvmatrix
      vnames <- paste("phi(", 1:p, ")", sep = "")
    }
    if(q > 0) {
      if(q > 1) {
        uumatrix <- (tccfAR(theta, theta)[ - (1:(q - 1))])[ - (
          q + 1)]
      }
      else if(q == 1) {
        uumatrix <- tccfAR(theta, theta)[ - (q + 1)]
      }
      uumatrix <- toeplitz(uumatrix)
      imatrix <- uumatrix
      unames <- paste("theta(", 1:q, ")", sep = "")
    }
    if(p > 0 && q > 0) {
      uvmatrix <- matrix(numeric(1), nrow = p, ncol = q)
      tuv <-  - tccfAR(phi, theta)
      for(i in 1:p) {
        for(j in 1:q) {
          uvmatrix[i, j] <- tuv[q + i - j]
        }
      }
      imatrix <- cbind(rbind(vvmatrix, t(uvmatrix)), rbind(uvmatrix, 
                                                           uumatrix))
    }
    inames <- c(vnames, unames)
    dimnames(imatrix) <- list(inames, inames)
    return(imatrix)
  }
