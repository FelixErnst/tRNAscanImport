#' @include tRNAscanImport.R
NULL

#' @name checktRNAscanGRanges
#' @aliases checktRNAscanGRanges
#' 
#' @title tRNAscan compatibility check
#' 
#' @description 
#' \code{checktRNAscanGRanges} checks whether a GRanges object contains the 
#' information expected for a tRNAscan result.
#' 
#' @param gr the \code{GRanges} object to test
#' 
#' @return a logical value
#' 
#' @examples 
#' file <- system.file("extdata", 
#'                     file = "sacCer3-tRNAs.ss.sort", 
#'                     package = "tRNAscanImport")
#' gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
#' checktRNAscanGRanges(gr)
NULL
#' @rdname checktRNAscanGRanges
#' @export
setMethod(
  f = "checktRNAscanGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr) .check_trna_granges(gr,
                                                TRNASCAN_FEATURES))

# checks whether a GRanges object is trnascan compatible
.check_trna_granges <- function(gr,features){
  if(class(gr) != "GRanges"){
    stop("Input is not a GRanges object.",
         call. = FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    stop("Input GRanges object does not meet the requirements of the ",
         "function. Please refer to the vignette of tRNAscanImport for ",
         "an exmaple on what information is expected.",
         call. = FALSE)
  }
  return(TRUE)
}
