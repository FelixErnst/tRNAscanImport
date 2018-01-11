#' @name tRNAScan2GRanges
#' 
#' @title tRNAScan2GRanges
#' 
#' \code{tRNAScan2GRanges} will convert a tRNAscan-SE output file into a 
#' GRanges object.
#'
#' @param file tRNAscan-SE input file
#' @param trim_intron optional logical: remove intron sequences
#'
#' @return a GRanges object
#' @export
#' 
#' @importFrom GenomicRanges GRanges
#' @importFrom assertive is_a_bool is_an_empty_string
#' @importFrom Biostrings DNAStringSet
#'
#' @examples
#' tRNAScan2GRanges(system.file("tRNAscan.txt",package="tRNAScan2GRanges"))
tRNAScan2GRanges <- function(file,
                              trim_intron = TRUE) {
  if(!assertive::is_a_bool(trim_intron)) trim_intron <- TRUE
  
  result <- .parse_tRNAscan(file)
  result <- lapply(result, function(trna){
    res <- list(no = trna$trna[3],
                chr = trna$trna[2])
    # If on minus strand
    if( as.numeric(trna$trna[5]) < as.numeric(trna$trna[4])){
      res <- append(res,
                    list(start = trna$trna[5],
                         end = trna$trna[4],
                         strand = "-"))  
    } else {
      res <- append(res,
                    list(start = trna$trna[4],
                         end = trna$trna[5],
                         strand = "+"))  
    }
    res <- append(res,
                  list(length = trna$trna[6],
                       type = trna$type[2],
                       ac = trna$type[3],
                       ac.start = trna$type[4],
                       ac.end = trna$type[5],
                       score = trna$type[6],
                       seq = trna$seq[2],
                       str = trna$str[2],
                       CCA.end = .has_CCA_end(trna$seq[2], trna$str[2]),
                       intron.start = trna$intron[4],
                       intron.end = trna$intron[5],
                       intron.locstart = trna$intron[2],
                       intron.locend = trna$intron[3],
                       hmm.score = trna$hmm[2],
                       sec.str.score = trna$hmm[3]))
    res <- lapply(res, function(x){if(is.null(x)) return(NA);x})
    return(res)
  })
  # create data.frame
  df <- setNames(lapply(names(result[[1]]), 
                        function(name){ unlist(lapply(result, 
                                                      function(x){unlist(x[[name]])})
                                               )}
                        ),
                 names(result[[1]]))
  df <- data.frame(df,
                  stringsAsFactors = FALSE)
  
  # optional: remove intron sequences
  if(trim_intron){
    df <- .cut_introns(df)
  }
  browser()
  # Contruct GRanges object
  gr <- GRanges(df)
  S4Vectors::mcols(gr)$seq <- DNAStringSet(S4Vectors::mcols(gr)$seq)
  
  return(gr)
}

# parse information on a tRNAscan file
.parse_tRNAscan <- function(file) {
  handle <- file(file, "r")
  res <- list()
  lines <- readLines(handle) 
  n <- 1
  repeat { 
    # read lines until empty line
    nextLines <- c()
    repeat {
      line <- lines[n]
      if(is.na(line) ||
         assertive::is_an_empty_string(line)) break
      nextLines <- c(nextLines,line)
      n <- n+1
    }
    if (length(nextLines) == 0 || n > length(lines)) break
    n <- n+1
    
    # parse the lines for tRNA data
    res <- append(res, list(.parse_tRNAscan_block(nextLines)))
  }
  close(handle)
  return(res)
}

# parse information on a single tRNA
.parse_tRNAscan_block <- function(lines) {
  .regex_custom <- function(line,regex_string){
    unlist(regmatches(line, 
                      regexec(regex_string, 
                              line)))
  }
  
  offset <- 0
  result <- list(trna = .regex_custom(lines[1], "([a-zA-Z0-9.:^*$@!+_?-|]+).trna([A-Z,-,_,0-9]+) \\(([0-9]+)-([0-9]+)\\).*Length: ([0-9]+) bp"),
                 type = .regex_custom(lines[2], "Type: ([A-z]{3}).*Anticodon: ([A-z]{3}) at ([0-9]+)-([0-9]+) .*Score: (.*)$"))
  
  intron <- .regex_custom(lines[3], "Possible intron: ([0-9]+)-([0-9]+) \\(([0-9]+)-([0-9]+)\\).*$")
  if(length(intron) != 0 ){
    offset <- offset+1
    result <- append(result, 
                     list(intron = intron))
  }
  hmm <- .regex_custom(lines[(3+offset)], "HMM Sc=([.,0-9]+).*Sec struct Sc=([.,0-9]+).*$")
  if(length(hmm) != 0 ){
    offset <- offset+1
    result <- append(result, 
                     list(hmm = hmm))
  }
  seq <- .regex_custom(lines[(4+offset)], "Seq: ([A,G,C,T,a,g,c,t]+)$")
  str <- .regex_custom(lines[(5+offset)], "Str: ([<,>,.]+)$")
  
  result <- append(result, list(seq = seq,
                                str = str))
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
  }
  return(FALSE)
}

# cuts out introns from sequence and structure
.cut_introns <- function(df){
  .cut_intron <- function(df, name){
    unlist(lapply(1:nrow(df), function(i){
      paste0(substring(df[i,name], 
                       1, 
                       as.numeric(tmp[i,]$intron.locstart)-1),
             substring(df[i,name], 
                       as.numeric(tmp[i,]$intron.locend)+1, 
                       nchar(tmp[i,name])))
    }))
  }
  
  tmp <- df[!is.na(df$intron.start),]
  seq <- .cut_intron(tmp,"seq")
  str <- .cut_intron(tmp,"str")
  
  df[!is.na(df$intron.start),"seq"] <- seq
  df[!is.na(df$intron.start),"str"] <- str
  return(df)
}

