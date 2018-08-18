#' @include tRNAscanImport.R
NULL

#' @rdname checktRNAscanGRanges
#' 
#' @export
setGeneric ( 
  name = "checktRNAscanGRanges",
  def = function(gr) standardGeneric("checktRNAscanGRanges")
)
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
  name = "gettRNAstructureGRanges",
  def = function(gr,
                 structure = "") standardGeneric("gettRNAstructureGRanges")
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