#' @include tRNA.R
NULL

STRUCTURE_OPEN_CHR <- c(">","\\[","\\(","\\{")
STRUCTURE_CLOSE_CHR <- c("<","\\]","\\)","\\}")

#' @name getBasePairing
#' @aliases getBasePairing gettRNABasePairing
#' 
#' @title Accessing Dot Bracket annotation
#' 
#' @description 
#' \code{getBasePairing} converts a dot bracket annotation into a
#' \code{data.frame}. Base pairing is indicated by corresponding numbers 
#' in the forward and reverse columns.
#' 
#' @return 
#' The result is a data.frame with three columns: pos, forward, reverse. If a
#' position is unpaired, forward and reverse will be \code{0}, otherwise it 
#' will match the base paired positions.
#' 
#' @param gr a GRanges object created by \code{import.tRNAscanAsGRanges} or 
#' GRanges with equivalent information. The \code{tRNA_str} column will be used
#' for input into \code{getBasePairing}.
#' 
#' @importFrom stringr str_locate_all
NULL

#' @rdname getBasePairing
#' 
#' @export
setMethod(
  f = "gettRNABasePairing",
  signature = signature(gr = "GRanges"),
  definition = function(gr) {
    .check_trnascan_granges(gr, TRNA_FEATURES)
    .get_base_pairing(gr$tRNA_str)
  }
)

#' @rdname getBasePairing
#' 
#' @param dotBracket character vectors describing a nucleotide sequence
#' structure in the dot bracket annotations. Valid characters are:
#' \code{.(\\{[><]\\})}
#' 
#' @export
getBasePairing <- function(dotBracket){
  assertive::assert_is_non_empty(dotBracket)
  .check_dot_bracket(dotBracket)
  dotBracket <- unlist(dotBracket)
  assertive::assert_all_are_non_missing_nor_empty_character(dotBracket)
  ans <- .get_base_pairing(dotBracket)
  return(ans)
}

# convert dot bracket annotation in ct like format
.get_base_pairing <- function(x){
  open <- lapply(STRUCTURE_OPEN_CHR,
                 function(chr){
                   stringr::str_locate_all(x, chr)
                 })
  close <- lapply(STRUCTURE_CLOSE_CHR,
                  function(chr){
                    stringr::str_locate_all(x, chr)
                  })
  lengthOpen <- lapply(open,function(z){lapply(z,length)})
  lengthClose <- lapply(close,function(z){lapply(z,length)})
  lengthMatch <- lapply(seq_along(lengthOpen),
                        function(i){
                          which(unlist(lengthOpen[[i]]) != unlist(lengthClose[[i]]))
                        })
  # check for unmatched positions
  if(any(unlist(lapply(lengthMatch,length)) != 0)){
    stop("Following structures are invalid: '",
         paste(unique(unlist(lengthMatch)),
               collapse = "'"),
         "' .\n They contain unmatched positions.")
  }
  structure <- mapply(.get_base_pairing_data_frame,
                      open,
                      close)
  structure <- split(structure,1:length(x))
  structure <- mapply(.complete_base_pairing_data_frame,
                      structure,
                      lapply(x,nchar),
                      SIMPLIFY = FALSE)
  return(structure)
}
# assembles base pairing data.frame
.get_base_pairing_data_frame <- function(open,
                                         close){
  lapply(seq_along(close),
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
}
# add missing value to base pairing data.frame
.complete_base_pairing_data_frame <- function(z,
                                              n){
  z <- do.call(rbind,z[!vapply(z,is.null,logical(1))])
  z2 <- z[,c(2,1)]
  colnames(z2) <- colnames(z)
  z <- rbind(z,z2)
  z$pos <- z$forward
  missing <- 1:n
  missing <- missing[!(missing %in% z$forward)]
  if(length(missing) > 0){
    z <- rbind(z,
               data.frame(pos = missing, 
                          forward = rep(0,length(missing)), 
                          reverse = rep(0,length(missing))))
  }
  z <- z[order(z$pos),c("pos","forward","reverse")]
  rownames(z) <- NULL
  return(z)
}