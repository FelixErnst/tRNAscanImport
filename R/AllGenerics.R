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

# Structures and Sequences -----------------------------------------------------

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
                 joinFeatures = FALSE,
                 padSequences = TRUE) standardGeneric("gettRNAstructureSeqs")
)
#' @rdname getBasePairing
#' 
#' @export
setGeneric ( 
  name = "gettRNABasePairing",
  def = function(gr) standardGeneric("gettRNABasePairing")
)

# Visualization ----------------------------------------------------------------

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