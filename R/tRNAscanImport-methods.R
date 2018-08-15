#' @include tRNAscanImport.R
NULL

#' @name gettRNAstructureSeq
#' 
#' @title Subsetting tRNA sequences
#'
#' @description
#' \code{gettRNAstructureSeq} returns split or partial tRNA sequences based on
#' the structure information of tRNAscan.
#' 
#' @param gr a GRanges object created by \code{import.tRNAscanAsGRanges} or 
#' GRanges with equivalent information.
#' @param structure optional parameter for returning just partial sequences.
#' The following values are accepted: 
#' \code{anticodonloop}, \code{Dloop}, \code{Tloop}, \code{acceptorStem}, 
#' \code{anticodonStem}, \code{DStem}, \code{TStem}, \code{variableLoop}, 
#' \code{discriminator}. (default = \code{structure = ""})
#' @param padSequences optional parameter whether sequences of the same type 
#' should be returned with the same length. For stems missing positions will be
#' filled up in the middle, for loops at the ends. 
#' (default = \code{padSequences = TRUE})
#' @param pad5prime If missing sequence length is an uneven number, 
#' should the 5'-end be padded by one more? (Only applies to loop sequences,
#' default = \code{padLeft = TRUE})
#'
#' @return 
#'  
#' 
#' @export
#' @examples
#' 
NULL

#' @rdname gettRNAstructureSeq
#' 
#' @export
setMethod(
  f = "gettRNAstructureSeq",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        structure,
                        padSequences,
                        pad5prime) {
    # input check
    .check_trnascan_granges(gr)
    .check_trna_structure_ident(structure)
    assertive::assert_is_a_bool(padSequences)
    assertive::assert_is_a_bool(pad5prime)
    if(structure == ""){
      structure <- TRNA_STRUCTURES
    }
    #
    
    browser()
    
    functionList <- list(anticodonloop = ".getAnticodonloop",
                         Dloop = ".getDloop",
                         Tloop = ".getTloop",
                         acceptorStem = ".getAcceptorStem",
                         anticodonStem = ".getAnticodonStem",
                         DStem = ".getDstem",
                         TStem = ".getTstem",
                         variableLoop = ".getVariableLoop",
                         discriminator = ".getDiscriminator")
    functions <- functionList[structure]
    res <- lapply(functions,
                  function(f,x){
                    pos <- do.call(f,x)
                    if(is.list(pos[[1]])){
                      return(lapply(pos,function(p){
                        IRanges::IRanges(p$start, p$end)
                      }))
                    }
                    if(is.list(pos)){
                      return(IRanges::IRanges(pos$start, pos$end))
                    }
                    return(IRanges(pos, pos))
                  },
                  list(gr))
    seqs <- lapply(res,
                   .assemble_sequences,
                   gr$tRNA_seq,
                   padSequences,
                   pad5prime)
    # gr[1]
    # switch(structure,
    #        FALSE = {
    #          
    #        },
    #        anticodonloop = ,
    #        Dloop = ,
    #        Tloop = ,
    #        acceptorStem = ,
    #        anticodonStem = ,
    #        DStem = ,
    #        TStem = ,
    #        variableLoop = ,
    #        discriminator = )
    # 
    if(length(seqs) == 1){
      seqs <- seqs[[1]]
    }
    return(seqs)
  }
)

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
      str_locate(substr(x[f]$tRNA_str, start + 2, .get_tRNA_length_wo_introns(x)),">")[,"start"]
    coord$start[f] <- start
    coord$end[f] <- end
  }
  return(coord)
}
.getAnticodonloop <- function(x){
  anticodon <- .getAnticodonStem(x)
  coord <- list(start = (anticodon$prime5$end + 1),
                end = (anticodon$prime3$start - 1))
  return(coord)
}
.getVariableLoop <- function(x){
  anticodon <- .getAnticodonStem(x)
  tstem <- .getTstem(x)
  coord <- list(start = (anticodon$prime3$end + 1),
                end = (tstem$prime5$start - 1))
  return(coord)
}
.getTloop <- function(x){
  tstem <- .getTstem(x)
  coord <- list(start = (tstem$prime5$end + 1),
                end = (tstem$prime3$start - 1))
  if(any(is.na(tstem$prime5$start))) {
    f <- which(is.na(tstem$prime5$start))
    acceptor <- .getAcceptorStem(x)
    end <- acceptor$prime3$start[f] - 1
    s <- reverse(substr(x[f]$tRNA_str, 1, end))
    # assumes minimum length of T of 3 nt
    start <- end + 2 - 
      str_locate(s,"<")[,"start"]
    coord$start[f] <- start
    coord$end[f] <- end
  }
  return(coord)
}
.getDiscriminator <- function(x){
  return(as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length_wo_introns(x)-3,
                           .get_tRNA_length_wo_introns(x))))
}


.getAcceptorStem <- function(x){
  end <- as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length_wo_introns(x)-3,
                           .get_tRNA_length_wo_introns(x)))
  coord <- list("prime5" = list(start = rep(1,length(x)),
                                end = rep(8,length(x))),
                "prime3" = list(start = (end-7),
                                end = (end-1)))
  return(coord)
}
.getDstem <- function(x){
  acceptor <- .getAcceptorStem(x)
  start5 <- acceptor$prime5$end + 1
  # D loop is expected to at least three nt long
  end5 <- start5 +
    str_locate(substr(x$tRNA_str, start5 + 2, .get_tRNA_length_wo_introns(x)),"\\.\\.\\.")[,"start"]
  start3 <- start5 + 1 +
    str_locate(substr(x$tRNA_str, start5 + 2, .get_tRNA_length_wo_introns(x)),"<")[,"start"]
  end3 <- start3 + (end5 - start5)
  coord <- list("prime5" = list(start = start5,
                                end = end5),
                "prime3" = list(start = start3,
                                end = end3))
  if(any(coord$prime5$start == (coord$prime5$end - 1))){
    f <- which(coord$prime5$start == (coord$prime5$end - 1))
    coord$prime5$start[f] <- NA
    coord$prime5$end[f] <- NA
    coord$prime3$start[f] <- NA
    coord$prime3$end[f] <- NA
  }
  return(coord)
}
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
    str_locate(substr(x$tRNA_str, start5 + 2, .get_tRNA_length_wo_introns(x)),"\\.\\.\\.\\.")[,"start"]
  start3 <- start5 + 1 +
    str_locate(substr(x$tRNA_str, start5 + 2, .get_tRNA_length_wo_introns(x)),"<")[,"start"]
  end3 <- start3 + (end5 - start5)
  coord <- list("prime5" = list(start = start5,
                                end = end5),
                "prime3" = list(start = start3,
                                end = end3))
  return(coord)
}
.getTstem <- function(x){
  acceptor <- .getAcceptorStem(x)
  end3 <- acceptor$prime3$start - 1
  s <- reverse(substr(x$tRNA_str, 1, end3))
  # T loop is expected to at least three nt long
  start3 <- end3 + 2 -
    str_locate(s,"\\.\\.\\.")[,"start"]
  end5 <- end3 + 1 -
    str_locate(s,">")[,"start"]
  start5 <- end5 - (end3 - start3)
  coord <- list("prime5" = list(start = start5,
                                end = end5),
                "prime3" = list(start = start3,
                                end = end3))
  if(any(coord$prime3$start == (coord$prime3$end - 1))){
    f <- which(coord$prime3$start == (coord$prime3$end - 1))
    coord$prime5$start[f] <- NA
    coord$prime5$end[f] <- NA
    coord$prime3$start[f] <- NA
    coord$prime3$end[f] <- NA
  }
  return(coord)
}

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

.assemble_sequences <- function(ir,
                                seqs,
                                padSequences,
                                pad5prime){
  # if it is a stem
  if(is.list(ir)){
    prime5 <- XVector::subseq(seqs, ir$prime5)
    prime3 <- XVector::subseq(seqs, ir$prime3)
    #
    if(!padSequences){
      return(Biostrings::xscat(
        prime5,
        prime3
      ))
    }
    #
    maxWidth <- max(width(prime5)+width(prime3))
    addN <- maxWidth - (width(prime5) + width(prime3))
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
    return(Biostrings::xscat(
      prime5,
      addString,
      prime3
    ))
  }
  # if it is a loop
  if(!padSequences){
    return(XVector::subseq(seqs, ir))
  }
  ans <- XVector::subseq(seqs, ir)
  maxWidth <- max(width(ans))
  missingWidth <- maxWidth - width(ans)
  addNLeft <- floor(missingWidth / 2)
  addNRight <- ceiling(missingWidth / 2)
  if(pad5prime){
    f <- addNRight > 0 & (missingWidth %% 2) == 1
    addNLeft[f] <- addNLeft[f] + 1
    addNRight[f] <- addNRight[f] - 1
  }
  addLeft <- DNAStringSet(unlist(
    lapply(addNLeft, 
           function(n){ 
             do.call(paste0,as.list(rep("-",(n))))
           })
  ))
  addRight <- DNAStringSet(unlist(
    lapply(addNRight, 
           function(n){ 
             do.call(paste0,as.list(rep("-",(n))))
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