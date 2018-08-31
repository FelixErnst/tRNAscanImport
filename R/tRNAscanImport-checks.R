#' @include tRNAscanImport.R
NULL

TRNASCAN_FEATURES <- c(
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

TRNA_STRUCTURES <- c(
  "anticodonloop",
  "Dloop",
  "Tloop",
  "acceptorStem",
  "anticodonStem",
  "DStem",
  "TStem",
  "variableLoop",
  "discriminator"
)

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
  definition = function(gr) .check_trnascan_granges(gr))

# checks whether a GRanges object is trnascan compatible
.check_trnascan_granges <- function(gr){
  if(class(gr) != "GRanges"){
    stop("Input is not a GRanges object.",
         call. = FALSE)
  }
  # check input
  if(length(intersect(TRNASCAN_FEATURES,colnames(S4Vectors::mcols(gr)))) !=
     length(TRNASCAN_FEATURES)){
    stop("Input GRanges object does not meet the requirements of the ",
         "function. Please refer to the vignette of tRNAscanImport for ",
         "an exmaple on what information is expected.",
         call. = FALSE)
  }
  return(TRUE)
}

# checks whether a string is a valid tRNA structure
.check_trna_structure_ident <- function(value,
                                        .xvalue = assertive::get_name_in_parent(value)){
  # check input
  checkValues <- c("",TRNA_STRUCTURES)
  if(!(value %in% checkValues)){
    stop("'",.xvalue,
         "' must be one of the following values: '",
         paste(checkValues, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}

# checks whether only dot bracket characters are present
.check_dot_bracket <- function(value,
                               .xvalue = assertive::get_name_in_parent(value)){
  checkChars <- c(STRUCTURE_OPEN_CHR,
                  STRUCTURE_CLOSE_CHR,
                  ".")
  checkChars <- gsub("\\\\","",checkChars)
  testChars <- unique(unlist(strsplit(value,"")))
  if(!all(testChars %in% checkChars)){
    prob <- testChars[!(testChars %in% checkChars)]
    stop("'",.xvalue,
         "' contains invalid characters for dot bracket annotation: '",
         paste(prob, collapse = "', '"),
         "'.",
         call. = FALSE)
  }
  return(invisible(TRUE))
}