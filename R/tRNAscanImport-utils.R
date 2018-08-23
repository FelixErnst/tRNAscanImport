#' @include tRNAscanImport.R
NULL

# get the tRNA length without the intron
.get_tRNA_length <- function(x){
  nchar(as.character(x$tRNA_seq))
}