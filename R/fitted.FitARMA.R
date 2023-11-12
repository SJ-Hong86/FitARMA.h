#' fitted method for class FitARMA
#' 
#' The fitted values are the observed minus residuals. If there is differencing, 
#' the observed values are those corresponding to the differenced time series.
#' 
#' @usage fitted(object, ...) ## S3 method for class 'FitARMA'
#' @param object class FitARMA object.
#' @param ... auxiliary parameters.
#' @returns vector or ts object.
#' @author A.I. McLeod.
#' @seealso [FitARMA()].
#' @examples 
#' data(SeriesA)
#' out<-FitARMA(SeriesA, c(1,0,1))
#' fitted(out)
#' 
#' @export
fitted.FitARMA <-
  function (object, ...) 
  {
    object$fits
  }
