#' print method for class FitARMA
#' 
#' a brief summary is printed out of the fitted model.
#' 
#' @usage print(x, ...) ## S3 method for class 'FitARMA'
#' @param x object, class FitARMA.
#' @param ... optional arguments.
#' @returns the result is displayed.
#' @author A.I. McLeod.
#' @seealso [FitARMA()].
#' @examples 
#' data(SeriesA)
#' out <- FitARMA(SeriesA, c(1,0,1))
#' print(out)
#' 
#' @export
print.FitARMA <-
  function (x, ...) 
  {
    LL <- x$loglikelihood
    k <- x$order[1]+x$order[3]
    if (!is.null(x$demean) && x$demean) 
      k <- k + 1
    n <- length(x$res)
    aic <- -2 * LL + 2 * k
    bic <- -2 * LL + log(n) * k
    dati <- x$DataTitle
    if (!is.null(dati)) 
      cat(dati, fill = TRUE)
    modti <- x$ModelTitle
    if (x$MeanMLE) 
      modti <- paste(modti, " With mean MLE.")
    cat(modti, fill = TRUE)
    cat(paste("length of series =", n, ",  number of parameters =", 
              k), fill = TRUE)
    cat(paste("loglikelihood =", round(LL, 2), ",  aic =", round(aic, 
                                                                 1), ",  bic = ", round(bic, 1)), fill = TRUE)
  }
