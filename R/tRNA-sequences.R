#' @include tRNA.R
NULL

#' @rdname gettRNAstructureSeqs
#'
#' @importFrom stringr str_locate
#' @importFrom BiocGenerics width
#' @importFrom IRanges reverse
#' @importFrom XVector subseq
#' @importFrom Biostrings DNAStringSet
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
    .check_trna_granges(gr, TRNA_FEATURES)
    .check_trna_structure_ident(structure)
    if(structure == ""){
      structure <- TRNA_STRUCTURES
    }
    assertive::assert_is_a_bool(joinCompletely)
    assertive::assert_is_a_bool(joinFeatures)
    assertive::assert_is_a_bool(padSequences)
    if(joinCompletely == TRUE && joinCompletely == joinFeatures){
      warning("Both 'joinCompletely' and 'joinFeatures' are set to TRUE.
              'joinCompletely' takes precedence.")
    }
    # Make sure sequences are a DNAStringSet
    if(class(gr$tRNA_seq) != "DNAStringSet"){
      gr$tRNA_seq <- Biostrings::DNAStringSet(gr$tRNA_seq)
    }
    # join completly or get splitup sequences
    if(joinCompletely){
      # get Ranges
      strList <- .get_base_pairing(gr$tRNA_str)
      res <- .get_tRNA_structures(tRNAStructureFunctionList,
                                  gr,
                                  strList)
      # get sequences
      seqs <- mapply(.assemble_sequences,
                     res,
                     names(res),
                     MoreArgs = list(gr,
                                     joinFeatures = FALSE,
                                     padSequences = TRUE,
                                     strList))
      # assemble boundaries IRanges
      z <- unlist(seqs)[TRNA_STRUCTURE_ORDER]
      ir <- lapply(lapply(unlist(seqs)[TRNA_STRUCTURE_ORDER],
                          BiocGenerics::width),
                   unique)
      start <- c(1,unlist(ir[1:(length(ir)-1)]))
      end <- unlist(ir)
      ir <- IRanges(start = unlist(lapply(seq_along(start),
                                          function(i){
                                            sum(start[1:i])
                                          })),
                    end = unlist(lapply(seq_along(end),
                                        function(i){
                                          sum(end[1:i])
                                        })))
      names(ir) <- names(end)
      # concat sequences
      seqs <- do.call(Biostrings::xscat,
                      z)
      # store boundaries as metadata
      S4Vectors::metadata(seqs) <- list("tRNA_structures" = ir)
    } else {
      # get Ranges
      strList <- gettRNABasePairing(gr)
      functions <- tRNAStructureFunctionList[structure]
      res <- .get_tRNA_structures(functions, gr, strList)
      seqs <- mapply(.assemble_sequences,
                     res,
                     names(res),
                     MoreArgs = list(gr,
                                     joinFeatures,
                                     padSequences,
                                     strList))
    }
    return(seqs)
  }
)

################################################################################
# join sequences of tRNA with correct padding
.assemble_sequences <- function(ir,
                                name,
                                gr,
                                joinFeatures,
                                padSequences,
                                strList){
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
                                  padSequences,
                                  strList))
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
    x <- .pad_unpaired_in_stem_region(seqs,ir,strList)
    prime5 <- .pad_right(x$prime5)
    prime3 <- .pad_left(x$prime3)
    ans <- list(prime5,
                prime3)
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if it is a stem and features should be joined and padded
  ##############################################################################
  if(is.list(ir) && joinFeatures){
    x <- .pad_unpaired_in_stem_region(seqs,ir,strList)
    prime5 <- .pad_right(x$prime5)
    prime3 <- .pad_left(x$prime3)
    maxWidth <- max(BiocGenerics::width(prime5) + BiocGenerics::width(prime3))
    addN <- maxWidth - (BiocGenerics::width(prime5) + BiocGenerics::width(prime3))
    addString <- Biostrings::DNAStringSet(rep("--",length(seqs)))
    addString <- Biostrings::xscat(
      addString,
      Biostrings::DNAStringSet(unlist(
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
    ans <- .pad_center(XVector::subseq(seqs, ir))
    names(ans) <- names(ir)
    return(ans)
  }
  ##############################################################################
  # if padding should be done left or right
  ##############################################################################
  if(name == "left" && padSequences){
    ans <- .pad_left(XVector::subseq(seqs, ir))
    names(ans) <- names(ir)
    return(ans)
  }
  if(name == "right" && padSequences){
    ans <- .pad_right(XVector::subseq(seqs, ir))
    names(ans) <- names(ir)
    return(ans)
  }
  #############################################
  # special rule for Dloop
  #############################################
  if(name == "Dloop" && padSequences){
    # three at 5'-end and one 3'-end
    # search for GGG, GG, GC, AG, AC, AT, TT, CT
    # if found put in der middle and the rest at the 3'-end
    # if not put everthing at the 5'-end
    GGfix <- c("GGG","GG","GC","AG","AC","AT","TT","CT")
    #
    getGGPos <- function(seqs, ir, i, searchString){
      if(is.null(searchString[i])) return(NULL)
      x <- stringr::str_locate(reverse(XVector::subseq(seqs, ir)),
                               reverse(searchString[i]))
      x <- BiocGenerics::start(ir) +
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
    #
    ir5 <- ir
    BiocGenerics::end(ir5) <- BiocGenerics::start(ir5) + 1
    ir3 <- ir
    BiocGenerics::start(ir3) <- BiocGenerics::end(ir3)
    irM <- ir
    BiocGenerics::start(irM) <- BiocGenerics::start(irM) + 2
    BiocGenerics::end(irM) <- BiocGenerics::end(irM) - 1
    #
    GGpos <- getGGPos(seqs, ir, 1, GGfix)
    BiocGenerics::start(irM[!is.na(GGpos)]) <- GGpos[!is.na(GGpos)]
    BiocGenerics::end(ir5[!is.na(GGpos)]) <- GGpos[!is.na(GGpos)] - 1
    # split in the middle if no GG like found
    if(any(is.na(GGpos))){
      BiocGenerics::start(irM[is.na(GGpos)]) <-
        BiocGenerics::start(irM[is.na(GGpos)]) +
        floor(BiocGenerics::width(irM[is.na(GGpos)])/2)
      BiocGenerics::end(ir5[is.na(GGpos)]) <-
        BiocGenerics::end(ir5[is.na(GGpos)]) +
        floor(BiocGenerics::width(irM[is.na(GGpos)])/2)
    }
    #
    ans <- Biostrings::xscat(
      .assemble_sequences(ir5,
                          "right",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE,
                          strList),
      .assemble_sequences(irM,
                          "right",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE,
                          strList),
      .assemble_sequences(ir3,
                          "left",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE,
                          strList)
    )
    #
    names(ans) <- names(ir)
    return(ans)
  }
  #############################################
  # special rule for variable loop
  #############################################
  if(name == "variableLoop" && padSequences){
    # keep one pos at 5'- and 3'-end
    ir5 <- ir
    BiocGenerics::end(ir5) <- BiocGenerics::start(ir5)
    prime5 <- XVector::subseq(seqs, ir5)
    ir3 <- ir
    BiocGenerics::start(ir3) <- BiocGenerics::end(ir3)
    prime3 <- XVector::subseq(seqs, ir3)
    # get the middle sequnce
    irM <- ir
    BiocGenerics::start(irM) <- BiocGenerics::start(irM) + 1
    BiocGenerics::end(irM) <- BiocGenerics::end(irM) - 1
    middle <- XVector::subseq(seqs, irM)
    # if longer variable loops exist search for paired regions
    facLong <- BiocGenerics::width(ir) > 4
    facBasePaired <- rep(FALSE, length(ir))
    if(any(facLong)){
      startPaired5 <- BiocGenerics::start(irM) - 1 +
        stringr::str_locate(substr(gr$tRNA_str,
                                   BiocGenerics::start(irM),
                                   BiocGenerics::end(irM)),
                            ">")[,"start"]
      facBasePaired <- !is.na(startPaired5)
      # if paired regions exist pad them accordingly
      if(any(facBasePaired)){
        startPaired5 <- BiocGenerics::start(irM[facBasePaired])
        # test for different loop types
        endPaired5 <- BiocGenerics::start(irM[facBasePaired]) - 1 +
          stringr::str_locate(substr(gr$tRNA_str[facBasePaired],
                                     BiocGenerics::start(irM[facBasePaired]),
                                     BiocGenerics::end(irM[facBasePaired])),
                              ">\\.\\.|>\\.<|><")[,"start"]
        startPaired3 <- BiocGenerics::start(irM[facBasePaired]) + 1 +
          stringr::str_locate(substr(gr$tRNA_str[facBasePaired],
                                     BiocGenerics::start(irM[facBasePaired]),
                                     BiocGenerics::end(irM[facBasePaired])),
                              "\\.\\.<|>><|\\.><")[,"start"]
        endPaired3 <- BiocGenerics::end(irM[facBasePaired])
        # assemble middle sequences
        stem <- list(prime5  = IRanges::IRanges(start = startPaired5,
                                                end = endPaired5),
                     prime3  = IRanges::IRanges(start = startPaired3,
                                                end = endPaired3))
        stem <- .assemble_sequences(stem,
                                    "variableloopstem",
                                    gr[facBasePaired],
                                    joinFeatures = FALSE,
                                    padSequences = TRUE,
                                    strList[facBasePaired])
        m <- XVector::subseq(seqs[facBasePaired], irM[facBasePaired])
        # proceed with sequences which have a loop
        f <- endPaired5 < startPaired3 - 1
        m[f] <- .assemble_sequences(IRanges::IRanges(start = endPaired5[f] + 1,
                                                     end = startPaired3[f] - 1),
                                    "right",
                                    gr[facBasePaired][f],
                                    joinFeatures = FALSE,
                                    padSequences = TRUE,
                                    strList[facBasePaired][f])
        # add spacer for missing loops
        addWidth <- rep(max(BiocGenerics::width(m[f])),length(m[!f]))
        m[!f] <- Biostrings::DNAStringSet(unlist(
          lapply(addWidth,
                 function(n){
                   do.call(paste0,as.list(rep("-",n)))
                 })
        ))
        # combien everything
        middle[facBasePaired] <- Biostrings::xscat(stem$prime5,m,stem$prime3)
      }
    }
    # pad non paired region sequences in the middle
    maxWidth <- max(BiocGenerics::width(middle))
    addNMiddle <- maxWidth - BiocGenerics::width(middle)
    addMiddle <- Biostrings::DNAStringSet(unlist(
      lapply(addNMiddle,
             function(n){
               do.call(paste0,as.list(rep("-",n)))
             })
    ))
    middle[addNMiddle > 0] <- Biostrings::xscat(
      XVector::subseq(middle[addNMiddle > 0],
                      1,
                      ceiling(width(middle[addNMiddle > 0])/2)),
      addMiddle,
      XVector::subseq(middle[addNMiddle > 0],
                      (ceiling(width(middle[addNMiddle > 0])/2) + 1),
                      width(middle[addNMiddle > 0]))
    )
    # assemble left right and middle sequences
    ans <- Biostrings::xscat(
      .assemble_sequences(ir5,
                          "right",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE,
                          strList),
      middle,
      .assemble_sequences(ir3,
                          "left",
                          gr,
                          joinFeatures = FALSE,
                          padSequences = TRUE,
                          strList)
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
  ans <- .pad_outside(XVector::subseq(seqs, ir))
  names(ans) <- names(ir)
  return(ans)
}

#
.pad_left <- function(ans){
  maxWidth <- max(BiocGenerics::width(ans))
  addNLeft <- maxWidth - BiocGenerics::width(ans)
  addLeft <- Biostrings::DNAStringSet(unlist(
    lapply(addNLeft,
           function(n){
             do.call(paste0,as.list(rep("-",n)))
           })
  ))
  ans[addNLeft > 0] <- Biostrings::xscat(
    addLeft,
    ans[addNLeft > 0]
  )
  return(ans)
}

#
.pad_right <- function(ans){
  maxWidth <- max(BiocGenerics::width(ans))
  addNRight <- maxWidth - BiocGenerics::width(ans)
  addRight <- Biostrings::DNAStringSet(unlist(
    lapply(addNRight,
           function(n){
             do.call(paste0,as.list(rep("-",n)))
           })
  ))
  ans[addNRight > 0] <- Biostrings::xscat(
    ans[addNRight > 0],
    addRight
  )
  return(ans)
}

#
.pad_center <- function(ans){
  maxWidth <- max(BiocGenerics::width(ans))
  addNMiddle <- maxWidth - BiocGenerics::width(ans)
  addMiddle <- Biostrings::DNAStringSet(unlist(
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
  return(ans)
}

#
.pad_outside <- function(ans){
  maxWidth <- max(BiocGenerics::width(ans))
  missingWidth <- maxWidth - BiocGenerics::width(ans)
  # if sequences should be padded not in the center but outside
  addNLeft <- floor(missingWidth / 2)
  addNRight <- ceiling(missingWidth / 2)
  f <- addNRight > 0 & (missingWidth %% 2) == 1
  addNLeft[f] <- addNLeft[f] + 1
  addNRight[f] <- addNRight[f] - 1
  addLeft <- Biostrings::DNAStringSet(unlist(
    lapply(addNLeft,
           function(n){
             do.call(paste0,as.list(rep("-",n)))
           })
  ))
  addRight <- Biostrings::DNAStringSet(unlist(
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
  return(ans)
}

# detect bulges in stem region and pad on opposite site
.pad_unpaired_in_stem_region <- function(seqs,
                                         ir,
                                         strList){
  prime5 <- XVector::subseq(seqs, ir$prime5)
  prime3 <- XVector::subseq(seqs, ir$prime3)
  if(length(unique(BiocGenerics::width(prime5))) > 1){
    f <- max(BiocGenerics::width(prime5)) > BiocGenerics::width(prime5)
    prime5[f] <- .add_padding_unpaired(prime5[f],
                                       ir$prime5[f],
                                       strList[f])
  }
  if(length(unique(BiocGenerics::width(prime3))) > 1){
    f <- max(BiocGenerics::width(prime3)) > BiocGenerics::width(prime3)
    prime3[f] <- .add_padding_unpaired(prime3[f],
                                       ir$prime3[f],
                                       strList[f])
  }
  return(list(prime5 = prime5,
              prime3 = prime3))
}

#
.add_padding_unpaired <- function(seqs,
                                  ir,
                                  strList){
  dims <- lapply(seq_along(strList),
                 function(i){
                   z <- strList[[i]][strList[[i]]$pos %in% BiocGenerics::start(ir[i]):BiocGenerics::end(ir[i]),]
                   z <- z[z$reverse != 0,]
                   zz <- vapply(seq_len(nrow(z)),
                                function(j){
                                  # missing pos on reverse
                                  z[j,]$reverse - 1 > z[j + 1,]$reverse &
                                    # not the last position
                                    j != nrow(z) &
                                    # not the same bulge on the other side
                                    (z[j,]$reverse - 1 - z[j + 1,]$reverse) != (z[j + 1,]$forward - 1 - z[j,]$forward)
                                },logical(1))
                   # if not unpaired position can be detected it is just a shorter
                   # stem
                   if(length(which(zz)) == 0){
                     return(NULL)
                   }
                   return(list(start = z[zz,]$forward + 1 - BiocGenerics::start(ir[i]),
                               stop = z[which(zz) + 1,]$forward +
                                 1 -
                                 BiocGenerics::start(ir[i]) -
                                 (z[which(zz) + 1,]$forward - z[which(zz),]$forward - 1),
                               length = z[zz,]$reverse -
                                 z[which(zz) + 1,]$reverse -
                                 1 -
                                 (z[which(zz) + 1,]$forward - z[which(zz),]$forward - 1)))
                 })
  dims <- dims[!vapply(dims,is.null,logical(1))]
  if(length(dims) == 0){
    return(seqs)
  }
  dims <- data.frame(dims)
  return(.add_padding_to_pos(seqs,
                             dims$start,
                             dims$stop,
                             dims$length))
}

#
.add_padding_to_pos <- function(seqs, start, stop, length){
  add <- Biostrings::DNAStringSet(unlist(
    lapply(length,
           function(n){
             do.call(paste0,as.list(rep("-",n)))
           })
  ))
  Biostrings::xscat(
    XVector::subseq(seqs, start = 1, end = start),
    add,
    XVector::subseq(seqs, start = stop, end = BiocGenerics::width(seqs))
  )
}
