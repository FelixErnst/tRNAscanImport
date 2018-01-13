#' @include AllGenerics-export.R
NULL

#' @name gettRNAscanSummary
#' @title tRNAscan2 Summary
#'
#' @aliases plottRNAscanSummary
#' 
#' @param gr a GRanges object created by the tRNAscan2GRanges or GRanges with
#' equivalent information. 
#'
#' @description
#' \code{gettRNAscanSummary()}: creates an DataFrame with aggregates of 
#' information
#' \code{plottRNAscanSummary()}: If ggplot2 is installed a plot of the data
#' is generated
#'
#' @return 
#' \code{gettRNAscanSummary()}: DataFrame
#' \code{plottRNAscanSummary()}: list of ggplots per column of data
#' 
#' @export
#' 
#' @importFrom S4Vectors DataFrame
#' @importFrom Biostrings alphabetFrequency
#' 
#' @examples
#' gr <- tRNAscan2GRanges(system.file("extdata", 
#'                              file = "sacCer3-tRNAs.ss.sort", 
#'                              package = "tRNAscan2GRanges"))
#' gettRNAscanSummary(gr)
#' plottRNAscanSummary(gr)
setMethod(
  f = "gettRNAscanSummary",
  signature = signature(gr = "GRanges"),
  definition = function(gr) {
    # check input
    checkCols <- c("length","type","anticodon","anticodon.start",
                   "anticodon.end","score","seq","str","CCA.end","intron.start",
                   "intron.end","intron.locstart","intron.locend","hmm.score",
                   "sec.str.score")
    if(length(intersect(checkCols,colnames(S4Vectors::mcols(gr)))) !=
       length(checkCols)){
      stop("Input GRanges object does not meet the requirements of the ",
           "function. Please refer to the vignette of tRNAscan2GRanges for ",
           "an exmaple on what information is expected.",
           call. = FALSE)
    }
    df <- S4Vectors::DataFrame(length = .get_lengths(gr),
                               gc = .get_gc_content(gr),
                               cca = .get_cca_ends(gr),
                               introns = .get_introns(gr),
                               .get_scores(gr))
    return(df)
  }
)

# check for equality of all elements of a list
.ident <- function(l){
  all(vapply(l, function(x){identical(x,l[[length(l)]])}, logical(1)))
}

# returns the length of tRNAs
.get_lengths <- function(gr){
  length <- as.numeric(S4Vectors::mcols(gr)$length)
  intron_locstart <- as.numeric(S4Vectors::mcols(gr)$intron.locstart)
  intron_locstart[is.na(intron_locstart)] <- 0
  intron_locend <- as.numeric(S4Vectors::mcols(gr)$intron.locend)
  intron_locend[is.na(intron_locend)] <- 0
  intron_length <- intron_locend - intron_locstart
  length <- length - vapply(intron_length, function(l){
    if(l > 0){
      return(l + 1)
    }
    0
  },numeric(1))
  length
}

# get GC content
.get_gc_content <- function(gr){
  freq <- Biostrings::alphabetFrequency(S4Vectors::mcols(gr)$seq, 
                                        baseOnly=TRUE)
  gc <- vapply(seq_len(nrow(freq)), function(i){
    sum(freq[i,c("G","C")]) / sum(freq[i,])
    },numeric(1))
  gc
}

# get frequence for the two first two bases in the amino acid stem of the tRNA
.get_freq_aa_stem_content <- function(gr){
  m <- Biostrings::oligonucleotideFrequency(
    Biostrings::subseq(S4Vectors::mcols(gr)$seq, 1, 2), 
    width = 2)
  df<-as.data.frame(m)
  aa <- reshape2::melt(df)
  colnames(aa) <- c("seq","value")
  aa
}

# fractions of tRNA with encoded CCA ends
.get_cca_ends <- function(gr){
  CCA <- S4Vectors::mcols(gr)$CCA.end
  as.numeric(CCA)
}

# fractions of tRNA with introns
.get_introns <- function(gr){
  introns <- S4Vectors::mcols(gr)$intron.start
  introns[is.na(introns)] <- 0
  introns[introns>0] <- 1
  introns
}

# aggregates the scores
.get_scores <- function(gr){
  as.list(S4Vectors::mcols(gr)[,c("score",
                                  "hmm.score",
                                  "sec.str.score")])
}


#' @rdname gettRNAscanSummary
#' 
#' @export
setMethod(
  f = "plottRNAscanSummary",
  signature = signature(grl = "GRangesList"),
  definition = function(grl) {
    if(!requireNamespace("ggplot2")){
      stop("ggplot2 is not installed.", call. = FALSE)
    }
    if(!requireNamespace("grDevices")){
      stop("ggplot2 is not installed.", call. = FALSE)
    }
    if(!requireNamespace("gridExtra")){
      stop("ggplot2 is not installed.", call. = FALSE)
    }
    
    if(length(grl) == 0)
      stop("GRangesList of length == 0 provided.",
           call. = FALSE)
    browser()
    data <- lapply(seq_len(length(grl)), function(i){
      
    })
    
    
    
    # don't plot during grid.arrange
    grDevices::pdf(file = NULL)
    
    plots <- lapply(data, function(plotData){
      plot <- .get_plot(plotData)
    })
    # don't plot during grid.arrange
    grDevices::dev.off()
    names(plots)
    return(plots)
  }
)

.get_plot <- function(df, name){
  browser()
  colNames <- colnames(df)
  df$id <- "tRNA"
  lapply(colNames, function(name){
    browser()
    if(name %in% c("cca","introns")){
      data <- df[,c("id",name)]
      colnames(data) <- c("id","value")
      data <- reshape2::dcast(data, id ~ value, fun.aggregate = length )
      data <- reshape2::melt(data, id = "id")
      return(ggplot2::ggplot(data, 
                             ggplot2::aes_(x = ~variable,
                                           y = ~value)) +
               ggplot2::geom_bar(stat = "identity") +
               ggplot2::coord_polar(theta="y"))
    }
    ggplot2::ggplot(df[,c("id",name)], 
                    ggplot2::aes_(x = ~id)) +
      ggplot2::geom_jitter(mapping = ggplot2::aes_string(y = name)) +
      ggplot2::scale_x_discrete(name = "") +
      ggplot2::ylim(0,NA)
  })
}


# 
# plotCCAEnds = function(tRNAs) {
#   lbls <- c("no CCA", "CCA")
#   percent <- signif(table(tRNAs$tRNA.CCA.end)/sum(table(tRNAs$tRNA.CCA.end))*100, digits=4)
#   lbls <- paste(lbls, percent) 
#   lbls <- paste(lbls,"%",sep="")
#   
#   pie(percent, labels=lbls, main="CCA ends encoded")
# }
# 
# plotGC = function(tRNAs){
#   GCs = table( substring( tRNAs$tRNA.seq, 0, 2 ) )
#   
#   lbls <- rownames(GCs)
#   percent <- round(GCs/sum(GCs)*100, digits=2)
#   lbls <- paste(lbls, percent) 
#   lbls <- paste(lbls,"%",sep="")
#   
#   pie(percent, labels=lbls, main="tRNA start sequence")
# }
