#' @include tRNAscanImport.R
NULL

#' @name istRNAscanGRanges
#' @aliases istRNAscanGRanges
#' 
#' @title tRNAscan compatibility check
#' 
#' @description 
#' \code{istRNAscanGRanges} checks whether a GRanges object contains the 
#' information expected for a tRNAscan result.
#' 
#' @param gr the \code{GRanges} object to test
#' 
#' @return a logical value
#' 
#' @examples 
#' file <- system.file("extdata", 
#'                     file = "yeast.tRNAscan", 
#'                     package = "tRNAscanImport")
#' gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
#' istRNAscanGRanges(gr)
NULL
#' @rdname istRNAscanGRanges
#' @export
setMethod(
  f = "istRNAscanGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr) .check_trnascan_granges(gr, TRNASCAN_FEATURES))

# checks whether a GRanges object is trnascan compatible
.check_trnascan_granges <- function(gr,features){
  if(!is(gr,"GRanges")){
    warning("Input is not a GRanges object.", call. = FALSE)
    return(FALSE)
  }
  # check input
  if(length(intersect(features,colnames(S4Vectors::mcols(gr)))) !=
     length(features)){
    warning("Input GRanges object does not meet the requirements of the ",
            "function. Please refer to the vignette of tRNAscanImport for ",
            "an exmaple on what information is expected.",
            call. = FALSE)
    return(FALSE)
  }
  return(TRUE)
}

#' @importClassesFrom BSgenome BSgenome
#' @importClassesFrom Rsamtools FaFile
.norm_genome <- function(genome){
  if(is(genome,"BSgenome")){
    return(genome)
  }
  if(is(genome,"DNAStringSet")){
    return(genome)
  }
  if(is.character(genome) && length(genome) == 1){
    genome <- try(Rsamtools::FaFile(genome))
  }
  if(is(genome,"FaFile")){
    Rsamtools::indexFa(genome)
    return(genome)
  }
  stop("'genome' must be an object of class 'BSgenome', 'DNAStringSet', 
       'FaFile' or character vector of length == 1L, which can be used to ",
       "create a 'FaFile' object.",
       call. = FALSE)
}
