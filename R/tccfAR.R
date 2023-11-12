#' Theoretical cross-covariances of auxilary AR process in ARMA(p,q)
#' 
#' The auxilary AR processes in the ARMA(p,q) model phi(B)z(t)=theta(B)a(t) are 
#' defined by phi(B)u(t)=-a(t) and theta(B)v(t)=a(t). The upper off-diagonal p-by-q 
#' block of the ARMA information matrix is obtained from the cross-covariances 
#' of u(t) and v(t). This function obtains these covariances. 
#' 
#' @usage tccfAR(phi, theta)
#' @param phi AR coefficients in ARMA.
#' @param theta MA coefficients in ARMA.
#' @details A set of linear equations which determine the covariances is solved. 
#'   The algorithm is similar in spirit to that for the autocovariances (McLeod, 1975).
#' @returns vector of cross-covariances.
#' @author A.I. McLeod.
#' @references McLeod, A.I. (1975), Derivation of the theoretical 
#'   autocorrelation function of autoregressive moving-average time 
#'   series, Applied Statistics 24, 255-256.
#' @seealso [InformationMatrixARMA()].
#' @examples
#' tccfAR(0.9,0.5)
#' 
#' @export
tccfAR <-
  function(phi, theta)
  {
    #auxilary function used with iarma#########
    #computes the theoretical cross-covariance function of two autoregressions
    # z[t]-phi[1] z_[t-1] --- phi[p] z[t-p]     = a[t]
    # z[t]-theta[1] z_[t-1] --- theta[q] z[t-q] = a[t]
    # where p, q are length(phi), length(theta)
    p <- length(phi)
    q <- length(theta)
    if(p == 0 || q == 0)
      return(numeric(0))
    k <- p + q
    rhs <- c(-1, rep(0, k - 1))
    A <- matrix(numeric(k^2), nrow = k, ncol = k)
    for(i in 1:k) {
      for(j in 1:k) {
        imj <- i - j
        ijq <- i + j - q - 1
        if(i > q) {
          if(i > j && imj <= q)
            A[i, j] <- theta[imj]
          else if(i > q && imj == 0)
            A[i, j] <- -1
        }
        else {
          if(ijq > 0 && ijq <= p)
            A[i, j] <- phi[ijq]
          else if(ijq == 0)
            A[i, j] <- -1
        }
      }
    }
    return(solve(A, rhs))
  }
