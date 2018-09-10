#' @include tRNAscanImport.R
NULL

################################################################################
# tRNAscanImport
################################################################################

#' @rdname istRNAscanGRanges
#' @export
setGeneric ( 
  name = "istRNAscanGRanges",
  def = function(gr) standardGeneric("istRNAscanGRanges")
)
#' @rdname gettRNAscanSummary
#' @export
setGeneric ( 
  name = "gettRNAscanSummary",
  def = function(gr) standardGeneric("gettRNAscanSummary")
) 

# Visualization ----------------------------------------------------------------

#' @rdname gettRNAscanSummary
#' @export
setGeneric ( 
  name = "plottRNAscan",
  def = function(grl) standardGeneric("plottRNAscan")
) 

#' @rdname gettRNAscanSummary
#' @export
setGeneric ( 
  name = "gettRNAscanPlots",
  def = function(grl) standardGeneric("gettRNAscanPlots")
) 


################################################################################
# tRNA
################################################################################

#' @rdname istRNAGRanges
#' @export
setGeneric ( 
  name = "istRNAGRanges",
  def = function(gr) standardGeneric("istRNAGRanges")
)

# Structures and Sequences -----------------------------------------------------

#' @rdname gettRNAstructureSeqs
#' @export
setGeneric (
  name = "gettRNAstructureGRanges",
  def = function(gr,
                 structure = "") standardGeneric("gettRNAstructureGRanges")
)
#' @rdname gettRNAstructureSeqs
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
#' @export
setGeneric (
  name = "gettRNABasePairing",
  def = function(gr) standardGeneric("gettRNABasePairing")
)

# Subsetting -------------------------------------------------------------------

#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasTStem",
  def = function(gr,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasTStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasDStem",
  def = function(gr,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasDStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasAcceptorStem",
  def = function(gr,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasAcceptorStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasAnticodonStem",
  def = function(gr,
                 length = NA,
                 unpaired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasAnticodonStem")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasTloop",
  def = function(gr,
                 length = NA) standardGeneric("hasTloop")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasDloop",
  def = function(gr,
                 length = NA) standardGeneric("hasDloop")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasAnticodonLoop",
  def = function(gr,
                 length = NA) standardGeneric("hasAnticodonLoop")
)
#' @rdname tRNA-subset
#' @export
setGeneric (
  name = "hasVariableLoop",
  def = function(gr,
                 length = NA,
                 paired = NA,
                 mismatches = NA,
                 bulged = NA) standardGeneric("hasVariableLoop")
)
