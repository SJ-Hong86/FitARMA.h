#' residuals method for class FitARMA
#' 
#' The innovation residuals are obtained.
#' 
#' @usage residuals(object, ...) ## S3 method for class 'FitARMA'
#' @param object class FitARMA object.
#' @param ... auxiliary parameters.
#' @returns vector or ts object.
#' @author A.I. McLeod.
#' @seealso [FitARMA()].
#' @examples 
#' data(SeriesA)
#' out <- FitARMA(SeriesA, c(1,0,1))
#' resid(out)
#' 
#' @export
residuals.FitARMA <-
  function (object, ...) 
  {
    object$res
  }
