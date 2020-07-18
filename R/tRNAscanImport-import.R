#' @include tRNAscanImport.R
NULL

#' @name import.tRNAscanAsGRanges
#' @aliases import.tRNAscanAsGRanges tRNAscan2GFF tRNAscanID
#' 
#' @title Importing a tRNAscan output file as a GRanges object
#' 
#' @description
#' The function \code{import.tRNAscanAsGRanges} will import a tRNAscan-SE output
#' file and return the information as a GRanges object. The reported 
#' intron sequences are spliced from the result by default, but can also 
#' returned as imported.
#' 
#' The function \code{tRNAScan2GFF} formats the output of 
#' \code{import.tRNAscanAsGRanges} to be GFF3 compliant.
#' 
#' \code{tRNAscanID} generates a unique tRNA ID, which is like the format used 
#' in the SGD annotation 
#' 
#' \code{t*AminoAcidSingleLetter*(*Anticodon*)*ChromosomeIdentifier**optionalNumberIfOnTheSameChromosome*}
#'  
#' Example: tP(UGG)L or tE(UUC)E1.
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
#' @param input 
#' \itemize{
#' \item \code{import.tRNAscanAsGRanges}: a tRNAscan-SE input file
#' \item \code{tRNAscan2GFF}: a compatible \code{GRanges} object such as the 
#' output of \code{import.tRNAscanAsGRanges}
#' }
#' @param as.GFF3 optional logical for \code{import.tRNAscanAsGRanges}: returns 
#' a gff3 compatible GRanges object directly. (default: \code{as.GFF3 = FALSE})
#' @param trim.intron optional logical for \code{import.tRNAscanAsGRanges}: 
#' remove intron sequences. This changes the tRNA length reported. To retrieve
#' the original length fo the tRNA gene, use the \code{width()} function on the 
#' GRanges object. (default: \code{trim.intron = TRUE})
#'
#' @return a GRanges object
#' @export
#' 
#' @importFrom Biostrings DNAStringSet
#' @importFrom Structstrings DotBracketStringSet
#' @importFrom BiocGenerics start
#' @importFrom GenomeInfoDb seqnames
#' @importFrom S4Vectors mcols
#' @importFrom stringr str_trim str_locate_all
#' @importFrom rtracklayer export.gff3
#'
#' @examples
#' gr <- import.tRNAscanAsGRanges(system.file("extdata", 
#'                                file = "yeast.tRNAscan", 
#'                                package = "tRNAscanImport"))
#' gff <- tRNAscan2GFF(gr)
#' identical(gff,import.tRNAscanAsGRanges(system.file("extdata", 
#'                                file = "yeast.tRNAscan", 
#'                                package = "tRNAscanImport"),
#'                                as.GFF3 = TRUE))
import.tRNAscanAsGRanges <- 
  function(input, as.GFF3 = FALSE, trim.intron = TRUE){
  # input check
  if(!.is_a_bool(as.GFF3)){
    warning("'as.GFF3' is not a bool. Resetting 'as.GFF3' == FALSE.",
            call. = FALSE)
    as.GFF3 <- FALSE
  }
  if(!.is_a_bool(trim.intron)){
    warning("'trim.intron' is not a bool. Resetting 'trim.intron' == TRUE.",
            call. = FALSE)
    trim.intron <- TRUE
  }
  # get tRNAscan as data.frame
  df <- .read_tRNAscan(input)
  # optional: remove intron sequences
  if(trim.intron){
    df <- .cut_introns(df)
  }
  # Contruct GRanges object
  gr <- GRanges(df)
  S4Vectors::mcols(gr)$tRNA_seq <- 
    Biostrings::DNAStringSet(S4Vectors::mcols(gr)$tRNA_seq)
  S4Vectors::mcols(gr)$tRNA_str <- 
    Structstrings::DotBracketStringSet(S4Vectors::mcols(gr)$tRNA_str)
  names(S4Vectors::mcols(gr)$tRNA_seq) <- as.character(seqnames(gr))
  names(S4Vectors::mcols(gr)$tRNA_str) <- as.character(seqnames(gr))
  S4Vectors::mcols(gr)$tRNA_length <- 
    nchar(as.character(S4Vectors::mcols(gr)$tRNA_seq))
  # sort GRanges object
  gr <- gr[order(GenomeInfoDb::seqnames(gr), BiocGenerics::start(gr))]
  # convert to gff3 compatible GRanges object
  if(as.GFF3){
    gr <- tRNAscan2GFF(gr)
  }
  return(gr)
}

# create data.frame from tRNAscan file
.read_tRNAscan <- function(file){
  # parse the information as a list of named lists
  result <- .parse_tRNAscan(file)
  # aggregate the data
  result <- lapply(result, 
                   function(trna){
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
                                        tRNA_anticodon.start = as.integer(trna$type[4]),
                                        tRNA_anticodon.end = as.integer(trna$type[5]),
                                        tRNAscan_score = as.numeric(trna$type[6]),
                                        tRNA_seq = as.character(trna$seq[2]),
                                        tRNA_str = as.character(trna$str[2]),
                                        tRNA_CCA.end = as.logical(.has_CCA_end(trna$seq[2], 
                                                                               trna$str[2])),
                                        # do not force type - optional data
                                        tRNAscan_potential.pseudogene = 
                                          ifelse(length(!is.na(trna$pseudogene[2])) != 0,
                                                 !is.na(trna$pseudogene[2]),
                                                 FALSE),
                                        tRNAscan_intron.start = trna$intron[4],
                                        tRNAscan_intron.end = trna$intron[5],
                                        tRNAscan_intron.locstart = trna$intron[2],
                                        tRNAscan_intron.locend = trna$intron[3],
                                        tRNAscan_hmm.score = trna$hmm[2],
                                        tRNAscan_sec.str.score = trna$secstruct[2],
                                        tRNAscan_infernal = trna$infernal[2]))
                     # if a field returns NULL because it is not set switch to NA, since this
                     # will persist for data.frame creation
                     res[vapply(res, is.null, logical(1))] <- NA
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
  # set data types
  rownames(df) <- NULL
  df$no <- as.integer(df$no)
  df$tRNA_length <- as.integer(df$tRNA_length)
  df$tRNAscan_potential.pseudogene <- 
    as.logical(df$tRNAscan_potential.pseudogene)
  df$tRNAscan_intron.start <- as.integer(df$tRNAscan_intron.start)
  df$tRNAscan_intron.end <- as.integer(df$tRNAscan_intron.end)
  df$tRNAscan_intron.locstart <- as.integer(df$tRNAscan_intron.locstart)
  df$tRNAscan_intron.locend <- as.integer(df$tRNAscan_intron.locend)
  df$tRNAscan_hmm.score <- as.numeric(df$tRNAscan_hmm.score)
  df$tRNAscan_sec.str.score <- as.numeric(df$tRNAscan_sec.str.score)
  df$tRNAscan_infernal <- as.numeric(df$tRNAscan_infernal)
  df[is.na(df$tRNA_type),"tRNA_type"] <- "Und"
  return(df)
}

# generates number for factor used for splitting
.get_factor_numbers <- function(x,i){
  if(is.na(x[(i+1)])) return(NULL)
  return(c(rep_len(i,(x[i + 1] - x[i])),.get_factor_numbers(x,(i + 1))))
}

# parse information on a tRNAscan file
.parse_tRNAscan <- function(file) {
  # open handle and read all lines
  handle <- file(file, "r")
  res <- list()
  lines <- readLines(handle) 
  on.exit(close(handle))
  if(length(lines) <= 1) stop("Empty file.", call. = FALSE)
  # determine empty line positions
  cuts <- unlist(lapply(seq_along(lines), function(i){
    if(stringr::str_trim(lines[i]) == "") return(i)
  }))
  # generate splitting factor
  if(is.null(cuts)){
    stop("Invalid format. No empty lines found as delimiter for individual",
         " tRNA entries. Please check the format of the input file.",
         call. = FALSE)
  }
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
  # valid tRNA block has 5 lines minimum
  if(length(lines) < 5) return(NULL)
  
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


#' @rdname import.tRNAscanAsGRanges
#'
#' @export
tRNAscan2GFF <- function(input) {
  .check_trnascan_granges(input, TRNASCAN_FEATURES)
  tRNAscan <- input
  # patch GRanges object with necessary columns for gff3 comptability
  S4Vectors::mcols(tRNAscan)$tRNA_seq <- 
    as.character(S4Vectors::mcols(tRNAscan)$tRNA_seq)
  S4Vectors::mcols(tRNAscan)$tRNA_str <- 
    as.character(S4Vectors::mcols(tRNAscan)$tRNA_str)
  # generate unique tRNA ID
  # SGD like format is used 
  # t*AminoAcidSingleLetter*(*Anticodon*)*ChromosomeIdentifier*
  # *optionalNumberIfOnTheSameChromosome*
  # Example: tP(UGG)L or tE(UUC)E1
  S4Vectors::mcols(tRNAscan)$ID <- tRNAscanID(tRNAscan)
  S4Vectors::mcols(tRNAscan)$type <- "tRNA"
  S4Vectors::mcols(tRNAscan)$type <- 
    as.factor(S4Vectors::mcols(tRNAscan)$type)
  S4Vectors::mcols(tRNAscan)$source <- "tRNAscan-SE"
  S4Vectors::mcols(tRNAscan)$source <- 
    as.factor(S4Vectors::mcols(tRNAscan)$source)
  S4Vectors::mcols(tRNAscan)$score <- NA
  S4Vectors::mcols(tRNAscan)$score <- 
    as.numeric(S4Vectors::mcols(tRNAscan)$score)
  S4Vectors::mcols(tRNAscan)$phase <- NA
  S4Vectors::mcols(tRNAscan)$phase <- 
    as.integer(S4Vectors::mcols(tRNAscan)$phase)
  S4Vectors::mcols(tRNAscan)$score <- 
    as.integer(S4Vectors::mcols(tRNAscan)$phase)
  # arrange columns in correct order
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

#' @rdname import.tRNAscanAsGRanges
#'
#' @export
tRNAscanID <- function(input){
  .check_trnascan_granges(input, TRNASCAN_FEATURES)
  tRNAscan <- input
  # create ids based on type, anticodon and chromosome
  chrom <- as.character(GenomeInfoDb::seqnames(tRNAscan))
  chromIndex <- unlist(lapply(seq_along(unique(chrom)), 
                              function(i){
                                rep(i,length(which(chrom == unique(chrom)[i])))
                              }))
  chromLetters <- .get_chrom_letters(length(unique(chromIndex)))
  # Modified version of AA code since the tRNAscan uses a slight deviation
  # from the one defined in Biostrings
  # swap values and names first
  aacode <- Biostrings::AMINO_ACID_CODE
  aacode_names <- names(aacode)
  aacode_values <- aacode
  aacode <- aacode_names
  names(aacode) <- aacode_values
  aacode <- c(aacode, 
              "SeC" = "U", 
              "Und" = "X",
              "fMe" = "M",
              "iMe" = "M",
              "Sup" = "X")
  aa <- unlist(lapply(tRNAscan$tRNA_type, 
                      function(type){
                        aacode[names(aacode) == type]
                      }))
  if( length(aa) != length(tRNAscan$tRNA_anticodon) ||
      length(aa) != length(chromLetters[chromIndex]) ){
    stop("Unknown tRNA type identifier: ",
         paste(tRNAscan$tRNA_type[!(tRNAscan$tRNA_type %in% names(aacode))],
               collapse = ", "),
         "\nKnown type identifier are: ",
         paste(names(aacode), collapse = "','"),
         "'. If this is a genuine identifier, please let us know.",
         call. = FALSE)
  }
  id <- paste0("t",
               aa,
               "(",
               tRNAscan$tRNA_anticodon,
               ")",
               chromLetters[chromIndex])
  # make ids unique if more than one tRNA of the same type is on the same 
  # chromosome
  uniqueID <- id[!duplicated(id)]
  uniqueID <- uniqueID[!(uniqueID %in% id[duplicated(id)])]
  pos <- match(unique(id[duplicated(id)]),id)
  dupID <- unlist(lapply(pos, function(i){
    x <- id[id == id[i]]
    ipos <- which(id == id[i])
    res <- paste0(x,seq(length(x)))
    names(res) <- ipos
    res
  }))
  id[as.numeric(names(dupID))] <- as.character(dupID)
  id
}

# get character values from "A" to "ZZZ" for example
.get_chrom_letters <- function(n){
  let <- list()
  add <- ""
  i <- 1
  while(n > 0){
    # get new letters and retrieve remaining n
    let[[i]] <- .get_chrom_letters2(n,add)
    n <- let[[i]]$remainder
    i <- i + 1
    # 
    batch <- (i-2) %/% length(LETTERS) + 1
    ln <- (i-1) - (batch-1)*length(LETTERS)
    if(n > 0){
      add <- let[[batch]][["let"]][ln]
    }
  }
  let <- unlist(lapply(seq_along(let), function(x){let[[x]][["let"]]}))
  let
}
.get_chrom_letters2 <- function(n, add = ""){
  ret <- list(let = c(),
              remaineder = 0)
  if(n > 0){
    l <- LETTERS[seq_len(n)]
    l <- l[!is.na(l)]
    let <- paste0(add,l)
    ret <- list(let = let,
                remainder = (n - length(let)))
  }
  ret
}
