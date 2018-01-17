#' @name tRNAscan2GRanges
#' 
#' @title tRNAscan2GRanges tRNAscan2GFF
#' 
#' @description
#' \code{tRNAScan2GRanges} will convert a tRNAscan-SE output file into a 
#' GRanges object. tRNAscan-SE 1.3.1 output is expected. Intron sequences are 
#' removed by default, but can also returned untouched.
#' \code{tRNAScan2GFF} formats the output of \code{tRNAScan2GRanges} to be GFF3
#' compliant
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
#' @importFrom rtracklayer export.gff3
#'
#' @examples
#' tRNAscan2GRanges(system.file("extdata", 
#'                              file = "sacCer3-tRNAs.ss.sort", 
#'                              package = "tRNAscan2GRanges"))
#'                              
#' tRNAscan2GFF(system.file("extdata", 
#'                          file = "sacCer3-tRNAs.ss.sort", 
#'                          package = "tRNAscan2GRanges"))
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
  S4Vectors::mcols(gr)$tRNA_seq <- DNAStringSet(S4Vectors::mcols(gr)$tRNA_seq)
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
                  list(tRNA_length = as.numeric(trna$trna[6]),
                       tRNA_type = as.character(trna$type[2]),
                       tRNA_anticodon = as.character(trna$type[3]),
                       tRNA_anticodon.start = as.numeric(trna$type[4]),
                       tRNA_anticodon.end = as.numeric(trna$type[5]),
                       tRNAscan_score = as.numeric(trna$type[6]),
                       tRNA_seq = as.character(trna$seq[2]),
                       tRNA_str = as.character(trna$str[2]),
                       tRNA_CCA.end = as.logical(.has_CCA_end(trna$seq[2], 
                                                         trna$str[2])),
                       # do not force type - optional data
                       tRNAscan_potential.pseudogene = !is.na(trna$pseudogene[2]),
                       tRNAscan_intron.start = trna$intron[4],
                       tRNAscan_intron.end = trna$intron[5],
                       tRNAscan_intron.locstart = trna$intron[2],
                       tRNAscan_intron.locend = trna$intron[3],
                       tRNAscan_hmm.score = trna$hmm[2],
                       tRNAscan_sec.str.score = trna$secstruct[2],
                       tRNAscan_infernal = trna$infernal[2]))
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
  df$tRNAscan_potential.pseudogene <- 
    as.logical(df$tRNAscan_potential.pseudogene)
  df$tRNAscan_intron.start <- as.numeric(df$tRNAscan_intron.start)
  df$tRNAscan_intron.end <- as.numeric(df$tRNAscan_intron.end)
  df$tRNAscan_intron.locstart <- as.numeric(df$tRNAscan_intron.locstart)
  df$tRNAscan_intron.locend <- as.numeric(df$tRNAscan_intron.locend)
  df$tRNAscan_hmm.score <- as.numeric(df$tRNAscan_hmm.score)
  df$tRNAscan_sec.str.score <- as.numeric(df$tRNAscan_sec.str.score)
  df$tRNAscan_infernal <- as.numeric(df$tRNAscan_infernal)
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

# regex wrapper
.regex_custom <- function(line,regex_string){
  unlist(regmatches(line, 
                    regexec(regex_string, 
                            line)))
}

# parse information on a single tRNA using regular expressions
.parse_tRNAscan_block <- function(lines) {
  # valid tRNA block has 6 lines minimum
  if(length(lines) <= 5) return(NULL)
  
  offset <- 0
  result <- list(trna = .regex_custom(lines[1], "([a-zA-Z0-9.:^*$@!+_?-|]+).trna([A-Z,-,_,0-9]+) \\(([0-9]+)-([0-9]+)\\).*Length: ([0-9]+) bp"),
                 type = .regex_custom(lines[2], "Type: ([A-z]{3}).*Anticodon: ([A-z]{3}) at ([0-9]+)-([0-9]+) .*Score: (.*)$"))
  intron <- .regex_custom(lines[3], "Possible intron: ([0-9]+)-([0-9]+) \\(([0-9]+)-([0-9]+)\\).*$")
  if(length(intron) != 0 ){
    offset <- offset+1
    result <- append(result, 
                     list(intron = intron))
  }
  pseudogene <- .regex_custom(lines[(3+offset)], "Possible (pseudogene): .*$")
  hmm <- .regex_custom(lines[(3+offset)], "HMM Sc=([-.,0-9]+).*$")
  secstruct <- .regex_custom(lines[(3+offset)], "Sec struct Sc=([-.,0-9]+).*$")
  infernal <- .regex_custom(lines[(3+offset)], "Infernal Sc=([-.,0-9]+).*$")
  if(length(pseudogene) != 0 | 
     length(hmm) != 0 | 
     length(secstruct) != 0 | 
     length(infernal) != 0 ){
    offset <- offset+1
    result <- append(result, 
                     list(pseudogene = pseudogene,
                          hmm = hmm,
                          secstruct = secstruct,
                          infernal = infernal))
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
                       (as.numeric(df[i,]$tRNAscan_intron.locstart)-1)),
             substring(df[i,name], 
                       (as.numeric(df[i,]$tRNAscan_intron.locend)+1), 
                       nchar(df[i,name])))
    }))
  }
  tmp <- df[!is.na(df$tRNAscan_intron.start),]
  seq <- .cut_intron(tmp,"tRNA_seq")
  str <- .cut_intron(tmp,"tRNA_str")
  df[!is.na(df$tRNAscan_intron.start),"tRNA_seq"] <- seq
  df[!is.na(df$tRNAscan_intron.start),"tRNA_str"] <- str
  return(df)
}


#' @rdname tRNAscan2GRanges
#'
#' @export
tRNAscan2GFF <- function(file,
                         trim_intron = TRUE) {
  tRNAscan <- tRNAscan2GRanges::tRNAscan2GRanges(file,trim_intron)
  S4Vectors::mcols(tRNAscan)$tRNA_seq <- 
    as.character(S4Vectors::mcols(tRNAscan)$tRNA_seq)
  S4Vectors::mcols(tRNAscan)$ID <- .create_tRNAscan_id(tRNAscan)
  S4Vectors::mcols(tRNAscan)$type <- "tRNA"
  S4Vectors::mcols(tRNAscan)$source <- "tRNAscan-SE"
  S4Vectors::mcols(tRNAscan)$score <- "."
  S4Vectors::mcols(tRNAscan)$phase <- "."
  S4Vectors::mcols(tRNAscan) <- 
    cbind(S4Vectors::mcols(tRNAscan)[,c("source",
                                        "type",
                                        "score",
                                        "phase",
                                        "ID")],
          S4Vectors::mcols(tRNAscan)[,-which(colnames(
            S4Vectors::mcols(tRNAscan)) %in% 
              c("source",
                "type",
                "score",
                "phase",
                "ID"))])
  return(tRNAscan)
}

.create_tRNAscan_id <- function(tRNAscan){
  browser()
  chrom <- as.character(GenomeInfoDb::seqnames(tRNAscan))
  chromIndex <- unlist(lapply(seq_along(unique(chrom)), 
                              function(i){
                                rep(i,length(which(chrom == unique(chrom)[i])))
                                }))
  chromLetters <- .get_chrom_letters(length(unique(chromIndex)))
  id <- paste0("t",
               tRNAscan$tRNA_type,
               "(",
               tRNAscan$tRNA_anticodon,
               ")",
               chrom)
  id
}

.get_chrom_letters <- function(n,add){
  let <- LETTERS[seq_len(n)]
  let <- let[!is.na(let)]
  if(n > 26){
    
  }
}