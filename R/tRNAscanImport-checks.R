#' @include tRNAscanImport.R
NULL

# checks whether a granges object is trnascan compatible
.check_trnascan_granges <- function(gr){
  # check input
  checkCols <- c(
    "tRNA_length",
    "tRNA_type",
    "tRNA_anticodon",
    "tRNA_anticodon.start",
    "tRNA_anticodon.end",
    "tRNAscan_score",
    "tRNA_seq",
    "tRNA_str",
    "tRNA_CCA.end",
    "tRNAscan_potential.pseudogene",
    "tRNAscan_intron.start",
    "tRNAscan_intron.end",
    "tRNAscan_intron.locstart",
    "tRNAscan_intron.locend",
    "tRNAscan_hmm.score",
    "tRNAscan_sec.str.score",
    "tRNAscan_infernal"
  )
  if(length(intersect(checkCols,colnames(S4Vectors::mcols(gr)))) !=
     length(checkCols)){
    stop("Input GRanges object does not meet the requirements of the ",
         "function. Please refer to the vignette of tRNAscan2GRanges for ",
         "an exmaple on what information is expected.",
         call. = FALSE)
  }
  
  
  
  
}

# check for equality of all elements of a list
.ident <- function(l){
  all(vapply(l, function(x){identical(x,l[[length(l)]])}, logical(1)))
}