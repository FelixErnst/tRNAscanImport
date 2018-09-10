#' @include tRNA.R
NULL

#' @name istRNAGRanges
#' @aliases istRNAGRanges
#'
#' @title tRNA compatibility check
#'
#' @description
#' \code{istRNAGRanges} checks whether a GRanges object contains the
#' information expected for a tRNA result.
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
#' istRNAGRanges(gr)
NULL
#' @rdname istRNAGRanges
#' @export
setMethod(
  f = "istRNAGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr) .check_trna_granges(gr,
                                                TRNA_FEATURES))

# checks whether a GRanges object is tRNA compatible
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
