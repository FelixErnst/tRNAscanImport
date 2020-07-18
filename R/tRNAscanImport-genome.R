#' @include tRNAscanImport.R
#' @include tRNAscanImport-checks.R
NULL

#' @name get.tRNAprecursor
#' @title Get tRNA precursor sequences
#' 
#' @description 
#' \code{get.tRNAprecursor} retrieves tRNA precursor sequences from genomic
#' sequences. The length of 5'- and 3'-overhangs can be specifid individually.
#' The output is checked for validity against the tRNA sequences as reported by
#' tRNAscan.
#' 
#' The chromosomes names of tRNAscan input and genome sequences must be 
#' compatible.
#'
#' @param input a compatible \code{GRanges} object
#' @param genome a \code{\link[BSgenome:BSgenome-class]{BSgenome}} object, a 
#' \code{\link[Biostrings:XStringSet-class]{DNAStringSet}} object, a 
#' \code{\link[Rsamtools:FaFile-class]{FaFile}} object or a character vector 
#' with a single value referncing a file, which can be coerced to a 
#' \code{FaFile} object.
#' @param add.5prime,add.3prime the length of overhangs as a single integer 
#' value each. (default \code{add.5prime = 50L})  
#' @param trim.intron \code{TRUE} or \code{FALSE}: Should intron sequences be
#' included in the precursor sequences? (default \code{trim.intron = FALSE})  
#'
#' @return a \code{DNAStringSet} object, containing the precursor sequences.
#' @export
#' 
#' @importFrom XVector subseq
#' @importFrom IRanges IRanges
#'
#' @examples
#' library(BSgenome.Scerevisiae.UCSC.sacCer3)
#' file <- system.file("extdata",
#'                     file = "yeast.tRNAscan",
#'                     package = "tRNAscanImport")
#' gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
#' genome <- getSeq(BSgenome.Scerevisiae.UCSC.sacCer3)
#' names(genome) <- c(names(genome)[-17],"chrmt")
#' get.tRNAprecursor(gr, genome)
#' # this produces an error since the seqnames do not match
#' \dontrun{
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' names(genome) <- c(names(genome)[-17],"chrmt")
#' get.tRNAprecursor(gr, genome)
#' }
#' # ... but it can also be fixed
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' seqnames(genome) <- c(seqnames(genome)[-17],"chrmt")
#' get.tRNAprecursor(gr, genome)
get.tRNAprecursor <- function(input, genome, add.5prime = 50L,
                              add.3prime = add.5prime, trim.intron = FALSE) {
  # input checks
  .check_trnascan_granges(input, TRNASCAN_FEATURES)
  genome <- .norm_genome(genome)
  if(!is.integer(add.5prime) || !is.integer(add.3prime)){
    stop("'add.5prime' and 'add.3prime' must be integer values.",call. = FALSE)
  }
  if(any(add.5prime < 0L) || any(add.3prime < 0L)){
    stop("All 'add.5prime' and 'add.3prime' values must >= 0L.",call. = FALSE)
  }
  if(length(add.5prime) > 1L){
    if(length(add.5prime) != length(input)){
      stop("'add.5prime' must be a single integer value or of the same length ",
           "as input.", call. = FALSE)
    }
  } else {
    add.5prime <- rep(add.5prime,length(input))
  }
  if(length(add.3prime) > 1L){
    if(length(add.3prime) != length(input)){
      stop("'add.3prime' must be a single integer value or of the same length ",
           "as input.", call. = FALSE)
    }
  } else {
    add.3prime <- rep(add.3prime,length(input))
  }
  if(!.is_a_bool(trim.intron)) {
    warning("'trim.intron' is not a bool. Resetting 'trim.intron' == TRUE.",
            call. = FALSE)
    trim.intron <- TRUE
  }
  # pre-check that input and genome are compatible
  chrNamesTSC <- as.character(seqnames(input))
  chrNamesGenome <- seqnames(seqinfo(genome))
  if(any(!(chrNamesTSC %in% chrNamesGenome))){
    stop("Not all chromosomes referenced in the tRNAscan data are present in ",
         "the genome data.",
         call. = FALSE)
  }
  # get precursor sequences
  strandM <- as.logical(strand(input) == "-")
  ir <- ranges(input)
  start(ir[!strandM]) <- start(ir[!strandM]) - add.5prime[!strandM]
  end(ir[!strandM]) <- end(ir[!strandM]) + add.3prime[!strandM]
  end(ir[strandM]) <- end(ir[strandM]) + add.5prime[strandM]
  start(ir[strandM]) <- start(ir[strandM]) - add.3prime[strandM]
  if(is(genome,"DNAStringSet")){
    m <- match(seqnames(input),chrNamesGenome)
    ans <- XVector::subseq(genome[m],ir)
    ans[strandM] <- reverseComplement(ans[strandM])
  } else if(is(genome,"FaFile")) {
    ir <- GenomicRanges::GRanges(seqnames = seqnames(input), ranges = ir,
                                 strand = strand(input))
    ans <- getSeq(genome, param = ir)
  } else if(is(genome,"BSgenome")) {
    ans <- getSeq(genome, seqnames(input), start(ir), end(ir), 
                  strand = strand(input))
  } else {
    stop("")
  }
  # post process precursor sequences
  # remove intron sequences if desired
  intronsPresent <- c("tRNAscan_intron.locstart") %in% colnames(mcols(input)) &
    any(!is.na(mcols(input)[["tRNAscan_intron.locstart"]]))
  if(intronsPresent && trim.intron){
    ans <- .remove_introns(ans, input, add.5prime, add.3prime)
  }
  # post-check that input and genome are compatible
  # remove introns for check only, if trim.intron == FALSE
  comp <- ans
  if(intronsPresent && !trim.intron){
    comp <- .remove_introns(comp, input, add.5prime, add.3prime)
  }
  # clip 5' and 3' additions
  ircomp <- IRanges::IRanges(add.5prime + 1L, lengths(comp) - add.3prime)
  comp <- XVector::subseq(comp,ircomp)
  # check that reconstructed sequences match the sequences reported in tRNAscan
  if(any(input$tRNA_seq != comp)){
    stop("Sequences from tRNAscan data and genomic sequences do not match. ",
         "The most probable cause are differences in coordinates between ",
         "tRNAscan and genome. Make sure that the tRNAscan data matches the ",
         "genome release.",
         call. = FALSE)
  }
  if(c("ID") %in% colnames(mcols(input))){
    names(ans) <- paste0("pre_",mcols(input)[["ID"]])
  } else {
    names(ans) <- paste0("pre_",seqnames(input),".tRNA",mcols(input)[["no"]])
  }
  ans
}

.remove_introns <- function(seq, input, add.5prime, add.3prime){
  withIntrons <- !is.na(mcols(input)[["tRNAscan_intron.locstart"]])
  if(!any(withIntrons)){
    return(seq)
  }
  intronStart <- mcols(input)[,"tRNAscan_intron.locstart"]
  intronEnd <- mcols(input)[,"tRNAscan_intron.locend"]
  ircomp5 <- IRanges::IRanges(1L,
                              add.5prime[withIntrons] + 
                                intronStart[withIntrons] - 1L)
  ircomp3 <- IRanges::IRanges(add.5prime[withIntrons] + intronEnd[withIntrons] +
                                1L,
                              lengths(seq[withIntrons]))
  seq[withIntrons] <- Biostrings::xscat(XVector::subseq(seq[withIntrons],ircomp5),
                                        XVector::subseq(seq[withIntrons],ircomp3))
  seq
}
