#' @include tRNAscanImport.R
NULL

#' @name gettRNAstructureSeqs
#' @aliases gettRNAstructureSeqs gettRNAstructureGRanges
#' 
#' @title Subsetting tRNA sequences
#'
#' @description
#' \code{gettRNAstructureSeq} returns split or partial tRNA sequences based on
#' the structure information of tRNAscan. Variances in length of certain 
#' structure features can be padded.
#' 
#' Padding works upon certain rules, which can be controlled via 
#' \code{padCenter} and \code{pad5prime}. Both parameters only affect 
#' padding of loops. With \code{padCenter = TRUE} length is adjusted by adding
#' "-" characters in the middle, whereas \code{padCenter = FALSE} adds them at
#' each end. If \code{pad5prime = TRUE} more padding characters are added on the
#' 5'-end if there are uneven numbers of missing positions.
#' 
#' @param gr a GRanges object created by \code{import.tRNAscanAsGRanges} or 
#' GRanges with equivalent information.
#' @param structure optional parameter for returning just partial sequences.
#' The following values are accepted: 
#' \code{anticodonloop}, \code{Dloop}, \code{Tloop}, \code{acceptorStem}, 
#' \code{anticodonStem}, \code{DStem}, \code{TStem}, \code{variableLoop}, 
#' \code{discriminator}. (default = \code{structure = ""})
#' @param joinCompletely Should the sequence parts, which are to be returned, be
#' joined into one sequence? (default = \code{joinCompletely = TRUE}))
#' Setting this to TRUE excludes \code{joinFeatures} be set to TRUE as well. In
#' addition, \code{joinCompletely = TRUE} uses automatically all sequence
#' structures.
#' @param joinFeatures Should the sequence parts, which are to be returned and
#' are from the same structure type, be joined into one sequence?
#' (default = \code{joinCompletely = FALSE})) Setting this to TRUE excludes 
#' \code{joinCompletely} be set to TRUE as well. \code{joinCompletely} takes
#' precedence.
#' @param padSequences parameter whether sequences of the same type 
#' should be returned with the same length. For stems missing positions will be
#' filled up in the middle, for loops at the ends. 
#' (default = \code{padSequences = TRUE}). If \code{joinCompletely == TRUE} this
#' is set to TRUE automatically.
#'
#' @return a list of \code{GRanges} or \code{DNAStringSet} objects. In case
#' joinCompletly is set to TRUE a single \code{DNAStringSet} is returned.
#' 
#' @export
#' @examples
#' gr <- import.tRNAscanAsGRanges(system.file("extdata", 
#'                                file = "sacCer3-tRNAs.ss.sort", 
#'                                package = "tRNAscanImport"))
#' gettRNAstructureGRanges(gr, structure = "anticodonloop")
#' gettRNAstructureSeqs(gr, structure = "anticodonloop")
NULL

tRNAStructureFunctionList <- list(
  anticodonloop = ".getAnticodonloop",
  Dloop = ".getDloop",
  Tloop = ".getTloop",
  acceptorStem = ".getAcceptorStem",
  anticodonStem = ".getAnticodonStem",
  DStem = ".getDstem",
  TStem = ".getTstem",
  variableLoop = ".getVariableLoop",
  discriminator = ".getDiscriminator")

#' @rdname gettRNAstructureSeqs
#' 
#' @importFrom IRanges IRanges
#' 
#' @export
setMethod(
  f = "gettRNAstructureGRanges",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        structure) {
    # input check
    .check_trnascan_granges(gr)
    .check_trna_structure_ident(structure)
    if(structure == ""){
      structure <- TRNA_STRUCTURES
    }
    #
    functions <- tRNAStructureFunctionList[structure]
    res <- lapply(functions,
                  function(f,x){
                    pos <- do.call(f,x)
                    if(is.list(pos[[1]])){
                      return(lapply(pos,function(p){
                        IRanges::IRanges(start = p$start,
                                         end = p$end,
                                         names = x[[1]]$tRNA_anticodon)
                      }))
                    }
                    if(is.list(pos)){
                      return(IRanges::IRanges(start = pos$start,
                                              end = pos$end,
                                              names = x[[1]]$tRNA_anticodon))
                    }
                    return(IRanges::IRanges(start = pos,
                                            end = pos,
                                            names = x[[1]]$tRNA_anticodon))
                  },
                  list(gr))
    return(res)
  }
)
#' @rdname gettRNAstructureSeqs
#' 
#' @importFrom stringr str_locate
#' @importFrom BiocGenerics width
#' @importFrom IRanges reverse
#' @importFrom XVector subseq
#' 
#' @export
setMethod(
  f = "gettRNAstructureSeqs",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        structure,
                        joinCompletely,
                        joinFeatures,
                        padSequences) {
    # input check
    assertive::assert_is_a_bool(joinCompletely)
    assertive::assert_is_a_bool(joinFeatures)
    assertive::assert_is_a_bool(padSequences)
    pad5prime <- TRUE
    if(joinCompletely == TRUE && joinCompletely == joinFeatures){
      warning("Both 'joinCompletely' and 'joinFeatures' are set to TRUE. 
              'joinCompletely' takes precedence.")
    }
    # Make sure sequences are a DNAStringSet
    if(class(gr$tRNA_seq) != "DNAStringSet"){
      gr$tRNA_seq <- DNAStringSet(gr$tRNA_seq)
    }
    #
    if(joinCompletely){
      # get Ranges
      res <- gettRNAstructureGRanges(gr, "")
      seqs <- mapply(.assemble_sequences,
                     res,
                     names(res),
                     MoreArgs = list(gr,
                                     joinFeatures = FALSE,
                                     padSequences = TRUE))
      seqs <- Biostrings::xscat(
        seqs$acceptorStem$prime5,
        seqs$DStem$prime5,
        seqs$Dloop,
        seqs$DStem$prime3,
        seqs$anticodonStem$prime5,
        seqs$anticodonloop,
        seqs$anticodonStem$prime3,
        seqs$variableLoop,
        seqs$TStem$prime5,
        seqs$Tloop,
        seqs$TStem$prime3,
        seqs$acceptorStem$prime3,
        seqs$discriminator
      )
    } else {
      # get Ranges
      res <- gettRNAstructureGRanges(gr, structure)
      seqs <- mapply(.assemble_sequences,
                     res,
                     names(res),
                     MoreArgs = list(gr,
                                     joinFeatures,
                                     padSequences))
    }
    return(seqs)
  }
)

####################################################
.getAcceptorStem <- function(x){
  end <- as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length_wo_introns(x)-3,
                           .get_tRNA_length_wo_introns(x)))
  coord <- list("prime5" = list(start = rep(1,length(x)),
                                end = rep(8,length(x))),
                "prime3" = list(start = unname(end-7),
                                end = unname(end-1)))
  return(coord)
}
####################################################
.getDstem <- function(x){
  acceptor <- .getAcceptorStem(x)
  # DStem starts for this with the last unpaired residue. However, the
  # difference towards the acceptor stem is also group towards the DStem 
  start5acceptor <- acceptor$prime5$end + 1
  start5 <- acceptor$prime5$end - 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               acceptor$prime5$end,
                               .get_tRNA_length_wo_introns(x)),
                        "\\.>")[,"start"]
  # D loop is expected to at least three nt long
  end5 <- start5 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length_wo_introns(x)),
                        "\\.\\.\\.")[,"start"]
  start3 <- start5 + 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length_wo_introns(x)),
                        "<")[,"start"]
  # if the T-stem is not 4 nt long, but the T loop is longer or equal than 6
  f <- (end5 - start5) < 4 & (start3 - end5) > 6
  if(any(f)){
    end5[f] <- end5[f] + 1
    start3[f] <- start3[f] - 1
  }
  #
  end3 <- start3 + (end5 - start5)
  coord <- list("prime5" = list(start = unname(start5acceptor),
                                end = unname(end5)),
                "prime3" = list(start = unname(start3),
                                end = unname(end3)))
  if(any(coord$prime5$start == (coord$prime5$end - 1))){
    f <- which(coord$prime5$start == (coord$prime5$end - 1))
    coord$prime5$start[f] <- NA
    coord$prime5$end[f] <- NA
    coord$prime3$start[f] <- NA
    coord$prime3$end[f] <- NA
  }
  return(coord)
}

####################################################
.getTstem <- function(x){
  acceptor <- .getAcceptorStem(x)
  end3 <- acceptor$prime3$start - 1
  s <- IRanges::reverse(substr(x$tRNA_str, 1, end3))
  # T loop is expected to at least three nt long
  start3 <- end3 + 2 -
    stringr::str_locate(s,
                        "\\.\\.\\.")[,"start"]
  end5 <- end3 + 1 - 3 -
    stringr::str_locate(s,
                        "\\.\\.\\.>")[,"start"]
  start5 <- end5 - (end3 - start3)
  coord <- list("prime5" = list(start = unname(start5),
                                end = unname(end5)),
                "prime3" = list(start = unname(start3),
                                end = unname(end3)))
  if(any(coord$prime3$start == (coord$prime3$end - 1))){
    f <- which(coord$prime3$start == (coord$prime3$end - 1))
    coord$prime5$start[f] <- NA
    coord$prime5$end[f] <- NA
    coord$prime3$start[f] <- NA
    coord$prime3$end[f] <- NA
  }
  return(coord)
}
####################################################
.getAnticodonStem <- function(x){
  dstem <- .getDstem(x)
  start <- dstem$prime3$end
  if(any(is.na(dstem$prime5$start))){
    f <- which(is.na(dstem$prime5$start))
    dloop <- .getDloop(x)
    start[f] <- dloop$end[f]
  }
  start5 <- start + 1
  # anticodon loop must at least be 4 nt long
  end5 <- start5 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length_wo_introns(x)),
                        "\\.\\.\\.\\.\\.")[,"start"]
  start3 <- start5 + 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length_wo_introns(x)),
                        "<")[,"start"]
  end3 <- start3 + (end5 - start5)
  # adjust coordinated for tRNA with non canonical base pairing at end of 
  # acceptor stem
  f <- (start3 - 9) > end5
  if(any(f)){
    end5[f] <- end5[f] + 1
    start3[f] <- start3[f] - 1
  }
  #
  coord <- list("prime5" = list(start = unname(start5),
                                end = unname(end5)),
                "prime3" = list(start = unname(start3),
                                end = unname(end3)))
  return(coord)
}
####################################################
####################################################
####################################################
.getDloop <- function(x){
  dstem <- .getDstem(x)
  coord <- list(start = (dstem$prime5$end + 1),
                end = (dstem$prime3$start - 1))
  if(any(is.na(dstem$prime5$start))) {
    f <- which(is.na(dstem$prime5$start))
    acceptor <- .getAcceptorStem(x)
    start <- acceptor$prime5$end[f]
    # assumes minimum length of D of 3 nt
    end <- start + 1 + 
      stringr::str_locate(substr(x[f]$tRNA_str,
                                 start + 2,
                                 .get_tRNA_length_wo_introns(x)),
                          ">")[,"start"]
    coord$start[f] <- start
    coord$end[f] <- end
  }
  return(coord)
}
####################################################
.getAnticodonloop <- function(x){
  anticodon <- .getAnticodonStem(x)
  coord <- list(start = (anticodon$prime5$end + 1),
                end = (anticodon$prime3$start - 1))
  return(coord)
}
####################################################
.getVariableLoop <- function(x){
  anticodon <- .getAnticodonStem(x)
  tstem <- .getTstem(x)
  coord <- list(start = (anticodon$prime3$end + 1),
                end = (tstem$prime5$start - 1))
  return(coord)
}
####################################################
.getTloop <- function(x){
  tstem <- .getTstem(x)
  coord <- list(start = (tstem$prime5$end + 1),
                end = (tstem$prime3$start - 1))
  if(any(is.na(tstem$prime5$start))) {
    f <- which(is.na(tstem$prime5$start))
    acceptor <- .getAcceptorStem(x)
    end <- acceptor$prime3$start[f] - 1
    s <- IRanges::reverse(substr(x[f]$tRNA_str, 1, end))
    # assumes minimum length of T of 3 nt
    start <- end + 2 - 
      stringr::str_locate(s,
                          "<")[,"start"]
    coord$start[f] <- start
    coord$end[f] <- end
  }
  return(coord)
}
####################################################
.getDiscriminator <- function(x){
  return(as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length_wo_introns(x)-3,
                           .get_tRNA_length_wo_introns(x))))
}

# get the tRNA length without the intron
.get_tRNA_length_wo_introns <- function(x){
  ans <- as.integer(x$tRNA_length) -
    (as.integer(ifelse(!vapply(x$tRNAscan_intron.locend,
                              is.na,
                              logical(1)),
                       x$tRNAscan_intron.locend,
                       0)) - 
    as.integer(ifelse(!vapply(x$tRNAscan_intron.locstart,
                             is.na,
                             logical(1)),
                      x$tRNAscan_intron.locstart,
                      0)))
  ans[!is.na(x$tRNAscan_intron.locend)] <- 
    ans[!is.na(x$tRNAscan_intron.locend)] - 1
  ans
}

#
.assemble_sequences <- function(ir,
                                name,
                                gr,
                                joinFeatures,
                                padSequences){
  seqs <- gr$tRNA_seq
  ##############################################################################
  # if it is a stem and features should be not be joined nor padded
  ##############################################################################
  if(is.list(ir) && !joinFeatures && !padSequences){
    ans <- mapply(.assemble_sequences,
                  ir,
                  MoreArgs = list(name,
                                  gr,
                                  joinFeatures,
                                  padSequences))
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should be joined, but not padded
  ##############################################################################
  if(is.list(ir) && joinFeatures && !padSequences){
    prime5 <- XVector::subseq(seqs, ir$prime5)
    prime3 <- XVector::subseq(seqs, ir$prime3)
    ans <- Biostrings::xscat(
      prime5,
      prime3
    )
    names(ans) <- name
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should not be joined, but padded
  ##############################################################################
  if(is.list(ir) && !joinFeatures && padSequences){
    prime5 <- XVector::subseq(seqs, ir$prime5)
    maxWidth5 <- max(BiocGenerics::width(prime5))
    addN5 <- maxWidth5 - BiocGenerics::width(prime5)
    addString5 <- DNAStringSet(unlist(
      lapply(addN5, 
             function(n){ 
               do.call(paste0,as.list(rep("-",n)))
             })
    ))
    prime5[addN5 > 0] <- Biostrings::xscat(
      prime5[addN5 > 0],
      addString5
    )
    prime3 <- XVector::subseq(seqs, ir$prime3)
    maxWidth3 <- max(BiocGenerics::width(prime3))
    addN3 <- maxWidth3 - BiocGenerics::width(prime3)
    addString3 <- DNAStringSet(unlist(
      lapply(addN3, 
             function(n){ 
               do.call(paste0,as.list(rep("-",n)))
             })
    ))
    prime3[addN3 > 0] <- Biostrings::xscat(
      addString3,
      prime3[addN3 > 0]
    )
    ans <- list(prime5,
                prime3)
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should be joined and padded
  ##############################################################################
  if(is.list(ir) && joinFeatures){
    prime5 <- XVector::subseq(seqs, ir$prime5)
    prime3 <- XVector::subseq(seqs, ir$prime3)
    maxWidth <- max(BiocGenerics::width(prime5)+BiocGenerics::width(prime3))
    addN <- maxWidth - (BiocGenerics::width(prime5) + BiocGenerics::width(prime3))
    addString <- DNAStringSet(rep("--",length(seqs)))
    addString <- Biostrings::xscat(
      addString, 
      DNAStringSet(unlist(
        lapply(addN, 
               function(n){ 
                 do.call(paste0,as.list(rep("-",(n + 1))))
               })
      ))
    )
    ans <- Biostrings::xscat(
      prime5,
      addString,
      prime3
    )
    names(ans) <- names(ir[[1]])
    return(ans)
  }
  ##############################################################################
  # special cases below
  ##############################################################################
  ##############################################################################
  # if padding should be done in the center
  ##############################################################################
  if(name == "center" && padSequences){
    ans <- XVector::subseq(seqs, ir)
    maxWidth <- max(BiocGenerics::width(ans))
    addNMiddle <- maxWidth - BiocGenerics::width(ans)
    addMiddle <- DNAStringSet(unlist(
      lapply(addNMiddle, 
             function(n){ 
               do.call(paste0,as.list(rep("-",n)))
             })
    ))
    ans[addNMiddle > 0] <- Biostrings::xscat(
      XVector::subseq(ans[addNMiddle > 0], 
                      1, 
                      ceiling(width(ans[addNMiddle > 0])/2)),
      addMiddle,
      XVector::subseq(ans[addNMiddle > 0],
                      (ceiling(width(ans[addNMiddle > 0])/2) + 1),
                      width(ans[addNMiddle > 0]))
    )
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if padding should be done left or right
  ##############################################################################
  # left is the default
  # if(name == "left" && padSequences){
  # }
  if(name == "right" && padSequences){
    ans <- XVector::subseq(seqs, ir)
    maxWidth <- max(BiocGenerics::width(ans))
    addNRight <- maxWidth - BiocGenerics::width(ans)
    addRight <- DNAStringSet(unlist(
      lapply(addNRight, 
             function(n){ 
               do.call(paste0,as.list(rep("-",n)))
             })
    ))
    ans[addNRight > 0] <- Biostrings::xscat(
      ans[addNRight > 0],
      addRight
    )
    names(ans) <- names(ir)
    return(ans)
  }
  #############################################
  # special rule for Dloop
  #############################################
  if(name == "Dloop" && padSequences){
    # three at 5'-end and one 3'-end
    # search for GG, GC, AG, AC, AT, TT, CT
    # if found put in der middle and the rest at the 3'-end
    # if not put everthing at the 5'-end
    GGfix <- c("GG","GC","AG","AC","AT","TT","CT")
    #
    getGGPos <- function(seqs, ir, i, searchString){
      if(is.null(searchString[i])) return(NULL)
      x <- stringr::str_locate(reverse(XVector::subseq(seqs, ir)),
                               reverse(searchString[i]))
      x <- start(ir) +
        width(reverse(XVector::subseq(seqs, ir))) -
        x[,"end"]
      names <-  rep(searchString[i], length(ir))
      names[is.na(x)] <- ""
      if(any(is.na(x))){
        y <- getGGPos(seqs[is.na(x)],
                      ir[is.na(x)],
                      (i + 1),
                      searchString)
        if(!is.null(y)){
          names[is.na(x)] <- names(y)
          x[is.na(x)] <- y
        }
      }
      names(x) <- names
      return(x)
    }
    # irInt <- ir
    # start(irInt) <- start(irInt) + 2
    # end(irInt) <- end(irInt) - 1
    # ans <- XVector::subseq(seqs, irInt)
    # aln <- msa::msa(ans, method = "ClustalW", gapOpening = 1, gapExtension = 1)
    # ans <- DNAStringSet(aln)
    ir5 <- ir
    end(ir5) <- start(ir5) + 1
    #
    ir3 <- ir
    start(ir3) <- end(ir3)
    #
    irM <- ir
    start(irM) <- start(irM) + 2
    end(irM) <- end(irM) - 1
    #
    GGpos <- getGGPos(seqs, ir, 1, GGfix)
    start(irM[!is.na(GGpos)]) <- GGpos[!is.na(GGpos)]
    end(ir5[!is.na(GGpos)]) <- GGpos[!is.na(GGpos)] - 1
    # split in the middle if no GG like found
    if(any(is.na(GGpos))){
      start(irM[is.na(GGpos)]) <- start(irM[is.na(GGpos)]) + 
        floor(width(irM[is.na(GGpos)])/2)
      end(ir5[is.na(GGpos)]) <- end(ir5[is.na(GGpos)]) + 
        floor(width(irM[is.na(GGpos)])/2)
    }
    #
    ans <- Biostrings::xscat(
      .assemble_sequences(ir5,
                          "right",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE),
      .assemble_sequences(irM,
                          "right",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE),
      .assemble_sequences(ir3,
                          "left",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE)
    )
    #
    names(ans) <- names(ir)
    return(ans)
  }
  #############################################
  # special rule for variable loop
  #############################################
  if(name == "variableLoop" && padSequences){
    # if only 3 or 4 long
    # if longer than 4 check for base pairing.
    # if base pairing put at the ends
    # if no base pairing put in the middle
    ir5 <- ir
    end(ir5) <- start(ir5)
    prime5 <- XVector::subseq(seqs, ir5)
    ir3 <- ir
    start(ir3) <- end(ir3)
    prime3 <- XVector::subseq(seqs, ir3)
    irM <- ir
    start(irM) <- start(irM) + 1
    end(irM) <- end(irM) - 1
    middle <- XVector::subseq(seqs, irM)
    browser()

    facLong <- width(ir) > 4
    facBasePaired <- rep(FALSE, length(ir))
    if(any(facLong)){
      startPaired5 <- start(ir) - 1 + 
        stringr::str_locate(substr(gr$tRNA_str,
                                   start(ir),
                                   end(ir)),
                            ">")[,"start"]
      facBasePaired <- !is.na(startPaired5)
      if(any(facBasePaired)){
        startPaired5 <- startPaired5[facBasePaired]
        # assumes at least three unpaired positions in loop
        endPaired5 <- start(ir[facBasePaired]) - 1 + 
          stringr::str_locate(substr(gr$tRNA_str[facBasePaired],
                                     start(ir[facBasePaired]),
                                     end(ir[facBasePaired])),
                              ">\\.")[,"start"]
        startPaired3 <- start(ir[facBasePaired]) + 2 + 
          stringr::str_locate(substr(gr$tRNA_str[facBasePaired],
                                     start(ir[facBasePaired]),
                                     end(ir[facBasePaired])),
                              "\\.\\.\\.<")[,"start"]
        endPaired3 <- end(ir[facBasePaired]) + 1 - 
          stringr::str_locate(reverse(substr(gr$tRNA_str[facBasePaired],
                                             start(ir[facBasePaired]),
                                             end(ir[facBasePaired]))),
                              "<")[,"start"]
        middle[facBasePaired] <- Biostrings::xscat(
          .assemble_sequences(IRanges::IRanges(start = ,
                                               end = ),
                              "right",
                              gr[facBasePaired],
                              joinFeatures = FALSE,
                              padSequences = TRUE)
        )
      }
    }
    middle[!facBasePaired] <- .assemble_sequences(irM[!facBasePaired],
                                                  "center",
                                                  gr[!facBasePaired],
                                                  joinFeatures = FALSE,
                                                  padSequences = TRUE)
    ans <- Biostrings::xscat(
      .assemble_sequences(ir5,
                          "right",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE),
      middle,
      .assemble_sequences(ir3,
                          "left",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE)
    )
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a loop and should not be padded
  ##############################################################################
  if(!is.list(ir) && !padSequences){
    ans <- XVector::subseq(seqs, ir)
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a loop and should be padded
  ##############################################################################
  ans <- XVector::subseq(seqs, ir)
  maxWidth <- max(BiocGenerics::width(ans))
  missingWidth <- maxWidth - BiocGenerics::width(ans)
  # if sequences should be padded not in the center but outside
  addNLeft <- floor(missingWidth / 2)
  addNRight <- ceiling(missingWidth / 2)
  f <- addNRight > 0 & (missingWidth %% 2) == 1
  addNLeft[f] <- addNLeft[f] + 1
  addNRight[f] <- addNRight[f] - 1
  addLeft <- DNAStringSet(unlist(
    lapply(addNLeft, 
           function(n){ 
             do.call(paste0,as.list(rep("-",n)))
           })
  ))
  addRight <- DNAStringSet(unlist(
    lapply(addNRight, 
           function(n){ 
             do.call(paste0,as.list(rep("-",n)))
           })
  ))
  ans[addNLeft > 0] <- Biostrings::xscat(
    addLeft,
    ans[addNLeft > 0]
  )
  ans[addNRight > 0] <- Biostrings::xscat(
    ans[addNRight > 0],
    addRight
  )
  names(ans) <- names(ir)
  return(ans)
}