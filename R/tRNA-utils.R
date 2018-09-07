#' @include tRNA.R
NULL

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

# get the tRNA length without the intron
.get_tRNA_length_c <- function(x){
  nchar(as.character(x$tRNA_seq))
}
.get_tRNA_length <- compiler::cmpfun(.get_tRNA_length_c)

.is_continous_evenly_spaced_c <- function(n){
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
.is_continous_evenly_spaced <- compiler::cmpfun(.is_continous_evenly_spaced_c)

# subset a structure data.frame using one or two iRanges
.subset_structure_c <- function(ir,
                                str,
                                pairedOnly = TRUE){
  if(pairedOnly){
    if(length(ir) == 2){
      f <- str$forward >= BiocGenerics::start(ir)[1] & 
        str$forward <= BiocGenerics::end(ir)[1] &
        str$reverse >= BiocGenerics::start(ir)[2] & 
        str$reverse <= BiocGenerics::end(ir)[2]
    } else {
      f <- str$forward >= BiocGenerics::start(ir) & 
        str$forward <= BiocGenerics::end(ir)
    }
  } else {
    if(length(ir) == 2){
      f <- (str$pos >= BiocGenerics::start(ir)[1] & 
              str$pos <= BiocGenerics::end(ir)[1]) |
        (str$pos >= BiocGenerics::start(ir)[2] & 
           str$pos <= BiocGenerics::end(ir)[2])
    } else {
      f <- str$pos >= BiocGenerics::start(ir) & 
        str$pos <= BiocGenerics::end(ir)
    }
  }
  str[f,]
}
.subset_structure <- compiler::cmpfun(.subset_structure_c)