#' @include tRNAscanImport.R
NULL

#' @rdname gettRNAscanSummary
#' 
#' @export
setGeneric ( 
  name = "gettRNAscanSummary",
  def = function(gr) standardGeneric("gettRNAscanSummary")
) 
#' @rdname gettRNAstructureSeq
#' 
#' @export
setGeneric ( 
  name = "gettRNAstructureSeq",
  def = function(gr,
                 structure = "",
                 padSequences = TRUE,
                 pad5prime = TRUE) standardGeneric("gettRNAstructureSeq")
) 

#' @rdname gettRNAscanSummary
#' 
#' @export
setGeneric ( 
  name = "plottRNAscan",
  def = function(grl) standardGeneric("plottRNAscan")
) 

#' @rdname gettRNAscanSummary
#' 
#' @export
setGeneric ( 
  name = "gettRNAscanPlots",
  def = function(grl) standardGeneric("gettRNAscanPlots")
) 