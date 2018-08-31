#' @include tRNAscanImport.R
NULL

# get the tRNA length without the intron
.get_tRNA_length <- function(x){
  nchar(as.character(x$tRNA_seq))
}

.is_continous_evenly_spaced <- function(n){
  if(length(n) < 2) return(FALSE)
  n <- n[order(n)]
  n <- n - min(n)
  step <- n[2] - n[1]
  test <- seq(from = min(n), to = max(n), by = step)
  if(length(n) == length(test) &&
     all(as.character(n) == as.character(test))){
    return(TRUE)
  }
  return(FALSE)
}

# subset a structure data.frame using one or two iRanges
.subset_strucutre <- function(ir5,
                              ir3,
                              str){
  if(!is.na(ir3)){
    f <- str$forward >= start(ir5) & str$forward <= end(ir5) &
      str$reverse >= start(ir3) & str$reverse <= end(ir3)
  } else {
    f <- str$forward >= start(ir5) & str$forward <= end(ir5)
  }
  str[f,]
}