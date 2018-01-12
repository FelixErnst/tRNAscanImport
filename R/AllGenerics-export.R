#' @include tRNAScan2GRanges.R
NULL

#' @rdname gettRNAscanSummary
#' 
#' @export
setGeneric ( 
  name = "gettRNAscanSummary",
  def = function(gr){standardGeneric("gettRNAscanSummary")} 
) 

#' @rdname gettRNAscanSummary
#' 
#' @export
setGeneric ( 
  name = "plottRNAscanSummary",
  def = function(gr){standardGeneric("plottRNAscanSummary")} 
) 