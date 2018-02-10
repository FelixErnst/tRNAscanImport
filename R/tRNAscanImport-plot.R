#' @include tRNAscanImport.R
NULL

#' @name gettRNAscanSummary
#' 
#' @title tRNAscanImport summary functions
#'
#' @aliases gettRNAscanSummary plottRNAscan gettRNAscanPlots
#' 
#' @param gr a GRanges object created by \code{import.tRNAscanAsGRanges} or 
#' GRanges with equivalent information. 
#' @param grl a GRangesList object created with GRanges created by 
#' \code{import.tRNAscanAsGRanges} or GRanges with equivalent information. 
#'
#' @description
#' \code{gettRNAscanSummary()}: creates an DataFrame with aggregated information
#' from the tRNAscan information.
#' 
#' \code{gettRNAscanPlots()}: If ggplot2 is installed a plot of the data
#' is generated and returned.
#' 
#' \code{plottRNAscan()}: If ggplot2 is installed a plot of the data
#' is generated and ploted directly.
#'
#' @return 
#' \code{gettRNAscanSummary()}: DataFrame
#' 
#' \code{gettRNAscanPlots()}: list of ggplots per column of data
#' 
#' \code{plottRNAscan()}:  directly plots the output of gettRNAscanPlots
#' 
#' 
#' @import methods
#' @importFrom S4Vectors DataFrame
#' @importFrom Biostrings alphabetFrequency
#' @importFrom reshape2 melt
#' 
#' @export
#' @examples
#' library(GenomicRanges, quietly = TRUE)
#' sce <- import.tRNAscanAsGRanges(system.file("extdata",
#'                                file = "sacCer3-tRNAs.ss.sort",
#'                                package = "tRNAscanImport"))
#' eco <- import.tRNAscanAsGRanges(system.file("extdata",
#'                         file = "eschColi_K_12_MG1655-tRNAs.ss.sort",
#'                         package = "tRNAscanImport"))
#' gettRNAscanSummary(sce)
#' plots <- gettRNAscanPlots(GRangesList(Sce = sce,
#'                                       Eco = eco))
NULL

#' @rdname gettRNAscanSummary
#' 
#' @export
setMethod(
  f = "gettRNAscanSummary",
  signature = signature(gr = "GRanges"),
  definition = function(gr) {
    .check_trnascan_granges(gr)
    df <- S4Vectors::DataFrame(length = .get_lengths(gr),
                               gc = .get_gc_content(gr),
                               cca = .get_cca_ends(gr),
                               pseudogene = .get_potential_pseudogene(gr),
                               introns = .get_introns(gr),
                               .get_scores(gr))
    return(df)
  }
)

# returns the length of tRNAs
.get_lengths <- function(gr){
  length <- as.numeric(S4Vectors::mcols(gr)$tRNA_length)
  intron_locstart <- as.numeric(S4Vectors::mcols(gr)$tRNAscan_intron.locstart)
  intron_locstart[is.na(intron_locstart)] <- 0
  intron_locend <- as.numeric(S4Vectors::mcols(gr)$tRNAscan_intron.locend)
  intron_locend[is.na(intron_locend)] <- 0
  intron_length <- intron_locend - intron_locstart
  length <- length - ifelse(intron_length > 0, intron_length + 1, 0)
  length
}

# get GC content
.get_gc_content <- function(gr){
  freq <- Biostrings::alphabetFrequency(S4Vectors::mcols(gr)$tRNA_seq, 
                                        baseOnly=TRUE)
  gc <- vapply(seq_len(nrow(freq)), function(i){
    sum(freq[i,c("G","C")]) / sum(freq[i,])
    },numeric(1))
  gc
}

# get frequence for the two first two bases in the amino acid stem of the tRNA
.get_freq_aa_stem_content <- function(gr){
  # requireNamespace("reshape2")
  m <- Biostrings::oligonucleotideFrequency(
    Biostrings::subseq(S4Vectors::mcols(gr)$tRNA_seq, 1, 2), 
    width = 2)
  df<-as.data.frame(m)
  aa <- reshape2::melt(df)
  colnames(aa) <- c("seq","value")
  aa
}

# fractions of tRNA with encoded CCA ends
.get_cca_ends <- function(gr){
  as.numeric(S4Vectors::mcols(gr)$tRNA_CCA.end)
}

# fractions of tRNA with pseudogene
.get_potential_pseudogene <- function(gr){
  as.numeric(S4Vectors::mcols(gr)$tRNAscan_potential.pseudogene)
}

# fractions of tRNA with introns
.get_introns <- function(gr){
  introns <- S4Vectors::mcols(gr)$tRNAscan_intron.start
  introns[is.na(introns)] <- 0
  introns[introns>0] <- 1
  introns
}

# aggregates the scores
.get_scores <- function(gr){
  as.list(S4Vectors::mcols(gr)[,c("tRNAscan_score",
                                  "tRNAscan_hmm.score",
                                  "tRNAscan_sec.str.score",
                                  "tRNAscan_infernal")])
}

#' @rdname gettRNAscanSummary
#' 
#' @export
setMethod(
  f = "plottRNAscan",
  signature = signature(grl = "GRangesList"),
  definition = function(grl) {
    plots <- gettRNAscanPlots(grl)
    if(!requireNamespace("graphics")){
      stop("Package 'graphics' is not installed.", call. = FALSE)
    }
    lapply(plots,graphics::plot)
    return(invisible(TRUE))
  }
)

#' @rdname gettRNAscanSummary
#' 
#' @export
setMethod(
  f = "gettRNAscanPlots",
  signature = signature(grl = "GRangesList"),
  definition = function(grl) {
    if(!requireNamespace("ggplot2")){
      stop("Package 'ggplot2' is not installed.", call. = FALSE)
    }
    if(!requireNamespace("scales")){
      stop("Package 'scales' is not installed.", call. = FALSE)
    }
    if(length(grl) == 0)
      stop("GRangesList of length == 0 provided.",
           call. = FALSE)
    
    # aggregate data
    data <- lapply(seq_len(length(grl)), function(i){
      mcoldata <- gettRNAscanSummary(grl[[i]])
      name <- names(grl[i])
      coldata <- lapply(seq_len(ncol(mcoldata)), function(i){
        column <- mcoldata[,i]
        column <- column[!is.na(column)]
        if(length(column) == 0) return(NULL)
        data.frame(id = name,
                   value = column)
      })
      names(coldata) <- colnames(mcoldata)
      return(coldata)
    })
    dataNames <- unique(unlist(lapply(data, names)))
    data <- lapply(dataNames, function(name){
      do.call(rbind, lapply(data, "[[", name))
    })
    names(data) <- dataNames
    # plot data
    plots <- lapply(seq_len(length(data)), function(i){
      if(is.null(data[[i]])){
        return(NULL)
      }
      .get_plot(data[i])
    })
    names(plots) <- dataNames
    plots <- plots[!vapply(plots, is.null, logical(1))]
    return(plots)
  }
)

# get a plot for one data type
.get_plot <- function(df){
  writtenNames <- list(length = "Length [nt]",
                       gc = "GC content [%]",
                       cca = "genomically encoded 3'-CCA ends [%]",
                       pseudogene = "Potential pseudogenes [%]",
                       introns = "Introns [%]",
                       tRNAscan_score = "tRNAscan-SE score",
                       tRNAscan_hmm.score = "HMM score",
                       tRNAscan_sec.str.score = "Secondary structure score",
                       tRNAscan_infernal = "Infernal score")
  dataType <- list(length = NA,
                   gc = "percent",
                   cca = "yn",
                   pseudogene = "yn",
                   introns = "yn",
                   tRNAscan_score = NA,
                   tRNAscan_hmm.score = NA,
                   tRNAscan_sec.str.score = NA,
                   tRNAscan_infernal = NA)
  name <- names(df)
  if(is.na(dataType[[name]])){
    plot <- ggplot2::ggplot(df[[name]],
                            ggplot2::aes_(x = ~id,
                                          y = ~value,
                                          colour = ~id)) +
      ggplot2::geom_violin(scale = "width") +
      ggplot2::geom_jitter(width = 0.2) +
      ggplot2::scale_x_discrete(name = "Organism") +
      ggplot2::scale_y_continuous(name = writtenNames[[name]]) +
      ggplot2::scale_colour_brewer(name = "Organism", palette = "Set1")
  } 
  if(!is.na(dataType[[name]]) && dataType[[name]] == "percent"){
    plot <- ggplot2::ggplot(df[[name]],
                    ggplot2::aes_(x = ~id,
                                  y = ~value,
                                  colour = ~id)) +
      ggplot2::geom_violin(scale = "width") +
      ggplot2::geom_jitter(width = 0.2) +
      ggplot2::scale_x_discrete(name = "Organism") +
      ggplot2::scale_y_continuous(name = writtenNames[[name]],
                                  breaks = c(0,0.25,0.5,0.75,1),
                                  labels = scales::percent,
                                  limits = c(0,1)) +
      ggplot2::scale_colour_brewer(name = "Organism", palette = "Set1")
  }
  if(!is.na(dataType[[name]]) && dataType[[name]] == "yn"){
    df[[name]][df[[name]]$value == 1,"colour"] <- "green"
    df[[name]][df[[name]]$value == 0,"colour"] <- "red"
    df[[name]][df[[name]]$value == 1,"value"] <- "Yes"
    df[[name]][df[[name]]$value == 0,"value"] <- "No"
    plot <- ggplot2::ggplot(df[[name]],
                            ggplot2::aes_(x = ~id,
                                          y = ~((..count..)/sum(..count..)),
                                          fill = ~colour)) +
      ggplot2::geom_bar(position = "fill") +
      ggplot2::scale_x_discrete(name = "Organism") +
      ggplot2::scale_y_continuous(name = writtenNames[[name]],
                                  labels = scales::percent,
                                  limits = c(0,1)) +
      ggplot2::scale_fill_identity(name = "",
                                   guide = "legend",
                                   labels = c("Yes","No"))
  }
  return(plot)
}