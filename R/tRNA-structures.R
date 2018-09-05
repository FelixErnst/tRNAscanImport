#' @include tRNA.R
NULL

#' @name gettRNAstructureSeqs
#' @aliases gettRNAstructureSeqs gettRNAstructureGRanges
#' 
#' @title tRNA structures and sequences
#'
#' @description
#' \code{gettRNAstructureGRanges} returns a list of GRanges describing the
#' boundaries of tRNA structures as extracted from the dot bracket annotation.
#' The dot bracket annotation is parsed using \code{gettRNABasePairing}, which
#' internally uses \code{getBasePairing}.
#' 
#' \code{gettRNAstructureSeq} returns split or partial tRNA sequences based on
#' the structure information of tRNAscan. Variances in length of certain 
#' structure features can be padded. If sequences are joined by setting
#' \code{joinCompletely = FALSE}, the boundaries of the tRNA structure are
#' stored in the result as metadata. They can be accessesed as an IRanges 
#' object by using \code{metadata()[["tRNA_structures"]]}.
#' 
#' @param gr a GRanges object created by \code{import.tRNAscanAsGRanges} or 
#' GRanges with equivalent information.
#' @param structure optional parameter for returning just partial sequences.
#' The following values are accepted: 
#' \code{anticodonloop}, \code{Dloop}, \code{Tloop}, \code{acceptorStem}, 
#' \code{anticodonStem}, \code{DStem}, \code{TStem}, \code{variableLoop}, 
#' \code{discriminator}. (default: \code{structure = ""})
#' @param joinCompletely Should the sequence parts, which are to be returned, be
#' joined into one sequence? (default: \code{joinCompletely = FALSE}))
#' Setting this to TRUE excludes \code{joinFeatures} be set to TRUE as well. In
#' addition, \code{joinCompletely = TRUE} uses automatically all sequence
#' structures.
#' @param joinFeatures Should the sequence parts, which are to be returned and
#' are from the same structure type, be joined into one sequence?
#' (default: \code{joinCompletely = FALSE})) Setting this to TRUE excludes 
#' \code{joinCompletely} be set to TRUE as well. \code{joinCompletely} takes
#' precedence.
#' @param padSequences parameter whether sequences of the same type 
#' should be returned with the same length. For stems missing positions will be
#' filled up in the middle, for loops at the ends. 
#' (default: \code{padSequences = TRUE}). If \code{joinCompletely == TRUE} this
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
#' gettRNABasePairing(gr[1:10])
#' getBasePairing(gr[1:10]$tRNA_str)
NULL

tRNAStructureFunctionList <- list(
  anticodonloop = ".getAnticodonloop",
  Dloop = ".getDloop",
  Tloop = ".getTloop",
  acceptorStem = ".getAcceptorStem",
  anticodonStem = ".getAnticodonStem",
  DStem = ".getDstem",
  TStem = ".getTstem",
  variableloop = ".getVariableLoop",
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
    .check_trnascan_granges(gr, TRNA_FEATURES)
    .check_trna_structure_ident(structure)
    if(structure == ""){
      structure <- TRNA_STRUCTURES
    }
    #
    functions <- tRNAStructureFunctionList[structure]
    strList <- gettRNABasePairing(gr)
    res <- .get_tRNA_structures(functions, gr, strList)
    return(res)
  }
)

.get_tRNA_structures <- function(functions,
                                 gr,
                                 strList){
  lapply(functions,
         function(f,x){
           pos <- do.call(f,list(x = x, strList = strList))
           if(is.list(pos[[1]])){
             return(lapply(pos,function(p){
               IRanges::IRanges(start = p$start,
                                end = p$end,
                                names = x$tRNA_anticodon)
             }))
           }
           if(is.list(pos)){
             return(IRanges::IRanges(start = pos$start,
                                     end = pos$end,
                                     names = x$tRNA_anticodon))
           }
           return(IRanges::IRanges(start = pos,
                                   end = pos,
                                   names = x$tRNA_anticodon))
         },
         gr)
}

####################################################
.getAcceptorStem <- function(x, strList){
  # acceptor stem must be 7 nt long
  # remove discriminator
  end <- as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length(x) - 3,
                           .get_tRNA_length(x))) - 1
  prime3 <- mapply(function(z,zz){z[(zz-8):(zz),]}, 
                   strList, 
                   end, 
                   SIMPLIFY = FALSE)
  last5primePaired <- lapply(prime3,
                             function(z){
                               z[z$reverse < 10 & z$reverse != 0,]$reverse[1]
                             })
  first3primePaired <- mapply(function(z,zz){z[z$reverse == zz,]$forward},
                              prime3,
                              last5primePaired,
                              SIMPLIFY = FALSE)
  coord <- list("prime5" = list(start = rep(1,length(x)),
                                # one position more is attributed to the 
                                # acceptor stem
                                end = unname(unlist(last5primePaired) + 1)), 
                "prime3" = list(start = unname(unlist(first3primePaired)),
                                # discriminator base is remove
                                end = end))
  return(coord)
}
####################################################
.getDstem <- function(x, strList){
  acceptor <- .getAcceptorStem(x, strList)
  # DStem starts for this with the last unpaired residue. However, the
  # difference towards the acceptor stem is also group towards the DStem 
  start5acceptor <- acceptor$prime5$end + 1
  start5 <- acceptor$prime5$end - 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               acceptor$prime5$end,
                               .get_tRNA_length(x)),
                        "\\.>")[,"start"]
  # D loop is expected to at least three nt long
  end5 <- start5 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length(x)),
                        "\\.\\.\\.")[,"start"]
  start3 <- start5 + 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length(x)),
                        "<")[,"start"]
  # if the Dstem is not 4 nt long, but the D loop is longer or equal than 6
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
.getTstem <- function(x, strList){
  acceptor <- .getAcceptorStem(x, strList)
  end3 <- acceptor$prime3$start - 1
  s <- IRanges::reverse(substr(x$tRNA_str, 1, end3))
  # T loop is expected to at least three nt long
  start3 <- end3 + 2 -
    stringr::str_locate(s,
                        "\\.\\.\\.")[,"start"]
  end5 <- end3 + 1 - 3 -
    stringr::str_locate(s,
                        "\\.\\.\\.>")[,"start"]
  start5 <- end3 + 1 - 
    stringr::str_locate(s,
                        ">\\.\\.")[,"start"]
  # if the T-stem is not 5 nt long on the 5'-site and 
  # the missing position at the 3'-end is unpaired include the missing 5'-pos 
  # as well
  f <- (end5 - start5) < 4 & (end3 - start3) > 3
  if(any(f)){
    start5[f] <- start5[f] - 1
  }
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
.getAnticodonStem <- function(x, strList){
  dstem <- .getDstem(x, strList)
  start <- dstem$prime3$end
  if(any(is.na(dstem$prime5$start))){
    f <- which(is.na(dstem$prime5$start))
    dloop <- .getDloop(x, strList)
    start[f] <- dloop$end[f]
  }
  start5 <- start + 1
  # anticodon loop must at least be 4 nt long
  end5 <- start5 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length(x)),
                        "\\.\\.\\.\\.\\.")[,"start"]
  start3 <- start5 + 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               start5 + 2,
                               .get_tRNA_length(x)),
                        "<")[,"start"]
  end3 <- start3 + 1 +
    stringr::str_locate(substr(x$tRNA_str,
                               start3 + 2,
                               .get_tRNA_length(x)),
                        "<\\.\\.|<\\.>|<>")[,"start"]
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
.getDloop <- function(x, strList){
  dstem <- .getDstem(x, strList)
  coord <- list(start = (dstem$prime5$end + 1),
                end = (dstem$prime3$start - 1))
  if(any(is.na(dstem$prime5$start))) {
    f <- which(is.na(dstem$prime5$start))
    acceptor <- .getAcceptorStem(x, strList)
    start <- acceptor$prime5$end[f]
    # assumes minimum length of D of 3 nt
    end <- start + 1 + 
      stringr::str_locate(substr(x[f]$tRNA_str,
                                 start + 2,
                                 .get_tRNA_length(x)),
                          ">")[,"start"]
    coord$start[f] <- start
    coord$end[f] <- end
  }
  return(coord)
}
####################################################
.getAnticodonloop <- function(x, strList){
  anticodon <- .getAnticodonStem(x, strList)
  coord <- list(start = (anticodon$prime5$end + 1),
                end = (anticodon$prime3$start - 1))
  return(coord)
}
####################################################
.getVariableLoop <- function(x, strList){
  anticodon <- .getAnticodonStem(x, strList)
  tstem <- .getTstem(x, strList)
  coord <- list(start = (anticodon$prime3$end + 1),
                end = (tstem$prime5$start - 1))
  return(coord)
}
####################################################
.getTloop <- function(x, strList){
  tstem <- .getTstem(x, strList)
  coord <- list(start = (tstem$prime5$end + 1),
                end = (tstem$prime3$start - 1))
  if(any(is.na(tstem$prime5$start))) {
    f <- which(is.na(tstem$prime5$start))
    acceptor <- .getAcceptorStem(x, strList)
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
.getDiscriminator <- function(x, strList){
  return(as.integer(ifelse(x$tRNA_CCA.end,
                           .get_tRNA_length(x)-3,
                           .get_tRNA_length(x))))
}