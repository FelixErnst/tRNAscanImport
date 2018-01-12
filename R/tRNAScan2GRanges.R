#' @name tRNAscan2GRanges
#' 
#' @title tRNAscan2GRanges
#' 
#' @description
#' \code{tRNAScan2GRanges} will convert a tRNAscan-SE output file into a 
#' GRanges object. tRNAscan-SE 1.3.1 output is expected. Intron sequences are 
#' removed by default, but can also returned untouched.
#'
#' @references 
#' Chan, Patricia P., and Todd M. Lowe. 2016. “GtRNAdb 2.0: An Expanded Database
#' of Transfer Rna Genes Identified in Complete and Draft Genomes.” Nucleic
#' Acids Research 44 (D1): D184–9. doi:10.1093/nar/gkv1309.
#'
#' Lowe, T. M., and S. R. Eddy. 1997. “TRNAscan-Se: A Program for Improved
#' Detection of Transfer Rna Genes in Genomic Sequence.” Nucleic Acids Research
#' 25 (5): 955–64. 
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
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_trim
#'
#' @examples
#' tRNAscan2GRanges(system.file("extdata", 
#'                              file = "tRNAscan.sort", 
#'                              package = "tRNAscan2GRanges"))
tRNAscan2GRanges <- function(file,
                              trim_intron = TRUE) {
  if(!assertive::is_a_bool(trim_intron)) trim_intron <- TRUE
    # get tRNAscan as data.frame
  df <- .read_tRNAscan(file)
  # optional: remove intron sequences
  if(trim_intron){
    df <- .cut_introns(df)
  }
  # Contruct GRanges object
  gr <- GRanges(df)
  S4Vectors::mcols(gr)$seq <- DNAStringSet(S4Vectors::mcols(gr)$seq)
  # sort GRanges object
  gr <- gr[order(GenomeInfoDb::seqnames(gr), BiocGenerics::start(gr))]
  return(gr)
}

# create data.frame from tRNAscan file
.read_tRNAscan <- function(file){
  # parse the information as a list of named lists
  result <- .parse_tRNAscan(file)
  # aggregate the data
  result <- lapply(result, function(trna){
    res <- list(no = as.numeric(trna$trna[3]),
                chr = as.character(trna$trna[2]))
    # If on minus strand
    if( as.numeric(trna$trna[5]) < as.numeric(trna$trna[4])){
      res <- append(res,
                    list(start = as.numeric(trna$trna[5]),
                         end = as.numeric(trna$trna[4]),
                         strand = "-"))  
    } else {
      res <- append(res,
                    list(start = as.numeric(trna$trna[4]),
                         end = as.numeric(trna$trna[5]),
                         strand = "+"))  
    }
    res <- append(res,
                  list(length = as.numeric(trna$trna[6]),
                       type = as.character(trna$type[2]),
                       anticodon = as.character(trna$type[3]),
                       anticodon.start = as.numeric(trna$type[4]),
                       anticodon.end = as.numeric(trna$type[5]),
                       score = as.numeric(trna$type[6]),
                       seq = as.character(trna$seq[2]),
                       str = as.character(trna$str[2]),
                       CCA.end = as.logical(.has_CCA_end(trna$seq[2], 
                                                         trna$str[2])),
                       # do not force type - optional data
                       intron.start = trna$intron[4],
                       intron.end = trna$intron[5],
                       intron.locstart = trna$intron[2],
                       intron.locend = trna$intron[3],
                       hmm.score = trna$hmm[2],
                       sec.str.score = trna$hmm[3]))
    # if a field returns NULL because it is not set switch to NA, since this
    # will persist for data.frame creation
    res <- lapply(res, function(x){if(is.null(x)) return(NA);x})
    return(res)
  })
  # create data.frame
  df <- lapply(names(result[[1]]), 
               function(name){ 
                 unlist(lapply(result, 
                               function(x){
                                 unlist(x[[name]])
                               }))
               })
  names(df) <- names(result[[1]])
  df <- data.frame(df,
                   stringsAsFactors = FALSE)
  df$intron.start <- as.numeric(df$intron.start)
  df$intron.end <- as.numeric(df$intron.end)
  df$intron.locstart <- as.numeric(df$intron.locstart)
  df$intron.locend <- as.numeric(df$intron.locend)
  df$hmm.score <- as.numeric(df$hmm.score)
  df$sec.str.score <- as.numeric(df$sec.str.score)
  return(df)
}

# generates number for factor used for splitting
.get_factor_numbers <- function(x,i){
  if(is.na(x[(i+1)])) return(NULL)
  return(c(rep_len(i,(x[i+1]-x[i])),.get_factor_numbers(x,(i+1))))
}

# parse information on a tRNAscan file
.parse_tRNAscan <- function(file) {
  # open handle and read all lines
  handle <- file(file, "r")
  res <- list()
  lines <- readLines(handle) 
  close(handle)
  if(length(lines) <= 1) stop("Empty file.", call. = FALSE)
  # determine empty line positions
  cuts <- unlist(lapply(seq_along(lines), function(i){
    if(stringr::str_trim(lines[i]) == "") return(i)
  }))
  # generate splitting factor
  cuts <- c(1,cuts,(max(cuts)+1))
  f <- as.factor(.get_factor_numbers(cuts,1))
  # parse the blocks for tRNA data
  res <- lapply(split(lines, f),function(nextLines){
    .parse_tRNAscan_block(nextLines[stringr::str_trim(nextLines) != ""])
  })
  res <- res[!vapply(res,is.null,logical(1))]
  if(length(res) == 0 ) stop("No tRNA information detected. Please make sure, ",
                             "that each tRNA is delimited by an empty line.",
                             call. = FALSE)
  return(res)
}

# parse information on a single tRNA using regular expressions
.parse_tRNAscan_block <- function(lines) {
  # valid tRNA block has 6 lines minimum
  if(length(lines) <= 5) return(NULL)
  # regex wrapper
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
    unlist(lapply(seq_len(nrow(df)), function(i){
      paste0(substring(df[i,name], 
                       1, 
                       (as.numeric(df[i,]$intron.locstart)-1)),
             substring(df[i,name], 
                       (as.numeric(df[i,]$intron.locend)+1), 
                       nchar(df[i,name])))
    }))
  }
  tmp <- df[!is.na(df$intron.start),]
  seq <- .cut_intron(tmp,"seq")
  str <- .cut_intron(tmp,"str")
  df[!is.na(df$intron.start),"seq"] <- seq
  df[!is.na(df$intron.start),"str"] <- str
  return(df)
}