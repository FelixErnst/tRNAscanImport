#' @name tRNAScan2gffAndFa
#' 
#' @title tRNAScan2gffAndFa
#' 
#' \code{tRNAScan2gffAndFa} will convert a tRNAscan-SE output file into a 
#' GRanges object.
#'
#' @param file tRNAscan-SE input file
#'
#' @return a GRanges object
#' @export
#' 
#' @import GenomicRanges
#'
#' @examples
#' tRNAScan2gffAndFa(system.file("tRNAscan.txt",package="tRNAScan2gffAndFa"))
tRNAScan2gffAndFa <- function(file, cluster = TRUE) {
  if(!is.logical(cluster)) cluster <- TRUE
  
  handle <- file(file, "r")
  
  tRNA.No <- c()
  tRNA.start <- c()
  tRNA.end <- c()
  tRNA.length <- c()
  tRNA.type <- c()
  tRNA.ac <- c()
  tRNA.ac.start <- c()
  tRNA.ac.end <- c()
  tRNA.score <- c()
  tRNA.seq <- c()
  tRNA.str <- c()
  tRNA.CCA.end <- c()
  
  repeat { 
    # read the lines 6 lines per tRNA
    lines <- readLines(handle, n = 6L) 
    if (length(lines) == 0) break 
    
    # parse the lines for tRNA data
    result = .parse_tRNA(lines)
    
    no <- length(tRNA.No)+1
    tRNA.No[no] <- result[[1]][2]
    tRNA.start[no] <- result[[1]][3]
    tRNA.end[no] <- result[[1]][4]
    tRNA.length[no] <- result[[1]][5]
    tRNA.type[no] <- result[[2]][2]
    tRNA.ac[no] <- result[[2]][3]
    tRNA.ac.start[no] <- result[[2]][4]
    tRNA.ac.end[no] <- result[[2]][5]
    tRNA.score[no] <- result[[2]][6]
    tRNA.seq[no] <- result[[3]][2]
    tRNA.str[no] <- result[[4]][2]
    # check if the tRNA has a CCA end encoded
    tRNA.CCA.end[no] <- .has_CCA_end(tRNA.seq[no], tRNA.str[no])
  }
  close(handle)
  
  df <- data.frame("tRNA.No" = tRNA.No,
                  "tRNA.start" = tRNA.start,
                  "tRNA.end" = tRNA.end,
                  "tRNA.length" = tRNA.length,
                  "tRNA.type" = tRNA.type,
                  "tRNA.ac" = tRNA.ac,
                  "tRNA.ac.start" = tRNA.ac.start,
                  "tRNA.ac.end" = tRNA.ac.end,
                  "tRNA.score" = tRNA.score,
                  "tRNA.seq" = tRNA.seq,
                  "tRNA.str" = tRNA.str,
                  "tRNA.CCA.end" = tRNA.CCA.end,
                  stringsAsFactors = FALSE)
  # df[] <- lapply(df, function(x) type.convert(as.character(x)))
  return(df)
}

# parse information on a single tRNA
.parse_tRNA <- function(lines) {
  result <- list(regmatches(lines[1], regexec(".+.trna([A-Z,-,_,0-9]+) \\(([0-9]+)-([0-9]+)\\).*Length: ([0-9]+) bp", lines[1])),
                 regmatches(lines[2], regexec("Type: ([A-z]{3}).*Anticodon: ([A-z]{3}) at ([0-9]+)-([0-9]+) .*Score: (.*)$", lines[2])),
                 regmatches(lines[4], regexec("Seq: ([A,G,C,T,a,g,c,t]+)$", lines[4])),
                 regmatches(lines[5], regexec("Str: ([<,>,.]+)$", lines[5])))
  return(result)
}

# check if a tRNA has a CCA end encoded
.has_CCA_end = function(seq, str) {
  end <- nchar( seq )
  start <- end - 2
  # last three nucleotides must be CCA and it must be unpaired
  if( substring( seq, start, end ) == "CCA" && 
      substring( str, start, end ) == "..." ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}