#' print method for class FitARMA
#' 
#' a summary is printed out of the fitted model
#' 
#' @usage summary(object, ...) ## S3 method for class 'FitARMA'
#' @param object class FitARMA object.
#' @param ... optional arguments.
#' @returns the result is displayed.
#' @author A.I. McLeod.
#' @seealso [FitARMA()].
#' @examples 
#' data(SeriesA)
#' out <- FitARMA(SeriesA, c(1,0,1))
#' summary(out)
#' 
#' @export
summary.FitARMA <-
  function (object, ...) 
  {
    LL <- object$loglikelihood
    k <- object$order[1]+object$order[3]
    if (!is.null(object$demean) && object$demean) 
      k <- k + 1
    n <- length(object$res)
    aic <- -2 * LL + 2 * k
    bic <- -2 * LL + log(n) * k
    dati <- object$DataTitle
    if (!is.null(dati)) 
      cat(dati, fill = TRUE)
    modti <- object$ModelTitle
    if (object$MeanMLE) 
      modti <- paste(modti, " With mean MLE.")
    cat(modti, fill = TRUE)
    cat(paste("length of series =", n, ",  number of parameters =", 
              k), fill = TRUE)
    cat(paste("loglikelihood =", round(LL, 2), ",  aic =", 
              round(aic, 1), ",  bic = ", round(bic, 1)), fill = TRUE)
  }
