#' @include tRNAscanImport.R
NULL

#' @rdname gettRNAscanSummary
#' 
#' @export
setGeneric ( 
  name = "gettRNAscanSummary",
  def = function(gr) standardGeneric("gettRNAscanSummary")
) 
#' @rdname gettRNAstructureSeqs
#' 
#' @export
setGeneric ( 
  name = "gettRNAstructureRanges",
  def = function(gr,
                 structure = "") standardGeneric("gettRNAstructureRanges")
)
#' @rdname gettRNAstructureSeqs
#' 
#' @export
setGeneric ( 
  name = "gettRNAstructureSeqs",
  def = function(gr,
                 structure = "",
                 joinCompletely = FALSE,
                 joinFeatures = TRUE,
                 padSequences = TRUE,
                 padCenter = TRUE,
                 pad5prime = FALSE) standardGeneric("gettRNAstructureSeqs")
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