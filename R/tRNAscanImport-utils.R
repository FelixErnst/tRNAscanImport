#' @include tRNAscanImport.R
NULL

# get the tRNA length without the intron
.get_tRNA_length <- function(x){
  nchar(as.character(x$tRNA_seq))
}

################################################################################

#' @rdname gettRNAstructureSeqs
#' 
#' @importFrom stringr str_locate_all
#' 
#' @export
setMethod(
  f = "gettRNABasePairing",
  signature = signature(gr = "GRanges"),
  definition = function(gr) {
    .check_trnascan_granges(gr)
    .get_base_pairing(gr$tRNA_str)
  }
)


#' @rdname gettRNAstructureSeqs
#' 
#' @export
getBasePairing <- function(dotBracket){
  assertive::assert_is_non_empty(dotBracket)
  .check_dot_bracket(dotBracket)
  dotBracket <- unlist(dotBracket)
  if(length(dotBracket) == 1){
    assertive::assert_is_a_non_empty_string(dotBracket)
    return(.get_base_pairing_per_str(dotBracket))
  }
  assertive::assert_all_are_non_missing_nor_empty_character(dotBracket)
  ans <- .get_base_pairing(dotBracket)
  return(ans)
}

# convert dot bracket annotation in ct like format
.get_base_pairing <- function(x){
  strList <- lapply(x,.get_base_pairing_per_str)
  return(strList)
}
.get_base_pairing_per_str <- function(str){
  open <- lapply(STRUCTURE_OPEN_CHR,
                 function(chr){
                   stringr::str_locate_all(str, chr)[[1]]
                 })
  close <- lapply(STRUCTURE_CLOSE_CHR,
                  function(chr){
                    stringr::str_locate_all(str, chr)[[1]]
                  })
  if(any(unlist(lapply(open,length)) != unlist(lapply(close,length)))){
    stop("Structure is invalid: ", str, ".\n It contains unmatched positions.")
  }
  structure <- lapply(seq_along(close),
                      function(i){
                        op <- open[[i]][,"start"]
                        cl <- close[[i]][,"start"]
                        if(length(cl) > 0){
                          forward <- rep(0, length(cl))
                          for(j in seq_along(cl)){
                            forward[j] <- rev(op[op < cl[j]])[1]
                            op <- op[op != forward[j]]
                          }
                          ans <- data.frame(forward = forward,
                                            reverse = cl)
                          return(ans)
                        }
                        return(NULL)
                      })
  structure <- structure[!vapply(structure,is.null,logical(1))]
  structure <- do.call(rbind,structure)
  structure2 <- structure[,c(2,1)]
  colnames(structure2) <- colnames(structure)
  structure <- rbind(structure,structure2)
  structure$pos <- structure$forward
  missing <- 1:nchar(str)
  missing <- missing[!(missing %in% structure$forward)]
  structure <- rbind(structure,
                     data.frame(pos = missing, 
                                forward = 0, 
                                reverse = 0))
  structure <- structure[order(structure$pos),c("pos","forward","reverse")]
  rownames(structure) <- NULL
  return(structure)
}