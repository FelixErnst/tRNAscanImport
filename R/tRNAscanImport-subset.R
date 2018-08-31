#' @include tRNAscanImport.R
NULL

#' @name tRNA-subset
#' 
#' @title Subsetting tRNA
#' 
#' @description 
#' title
#' 
#' @param gr a GRanges object from a tRNAscan import or with equivalent
#' information
#' 
#' @return logical vector
#' 
NULL

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasTStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        unpaired,
                        length,
                        mismatches,
                        bulged) .subset_tRNA_stem("TStem",
                                                  gr,
                                                  unpaired,
                                                  length,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasDStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        unpaired,
                        length,
                        mismatches,
                        bulged) .subset_tRNA_stem("DStem",
                                                  gr,
                                                  unpaired,
                                                  length,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasAcceptorStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        unpaired,
                        length,
                        mismatches,
                        bulged) .subset_tRNA_stem("acceptorStem",
                                                  gr,
                                                  unpaired,
                                                  length,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasAnticodonStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        unpaired,
                        length,
                        mismatches,
                        bulged) .subset_tRNA_stem("anticodonStem",
                                                  gr,
                                                  unpaired,
                                                  length,
                                                  mismatches,
                                                  bulged))

.subset_tRNA_stem <- function(ident,
                              gr,
                              unpaired,
                              length,
                              mismatches,
                              bulged){
  ans <- rep(TRUE,length(gr))
  # input check
  .check_trnascan_granges(gr, TRNASCAN_FEATURES)
  if(!(ident %in% names(tRNAStructureFunctionList))){
    stop("Unknown identifier '",ident,"'.")
  }
  if(!is.na(mismatches)){
    assertive::assert_is_a_bool(mismatches)
    unpaired <- TRUE
  }
  if(!is.na(bulged)){
    assertive::assert_is_a_bool(bulged)
    unpaired <- TRUE
  }
  #
  strList <- gettRNABasePairing(gr)
  str <- .get_tRNA_structures(tRNAStructureFunctionList[[ident]],
                              gr,
                              strList)[[1]]
  if(!is.na(unpaired) ||
     !is.na(length)){
    strList <- mapply(.subset_strucutre,
                      split(str$prime5,1:length(str$prime5)),
                      split(str$prime3,1:length(str$prime3)),
                      strList,
                      SIMPLIFY = FALSE)
  }
  browser()
  if(!is.na(unpaired)){
    assertive::assert_is_a_bool(unpaired)
    if(is.na(mismatches) && is.na(bulged)){
      pairedAns <- vapply(strList,
                          function(str){
                            .is_continous_evenly_spaced(str[str$forward != 0]$forward) &
                              .is_continous_evenly_spaced(str[str$reverse != 0]$reverse)
                          }
                          logical(1))
      if(assertive::is_true(unpaired)){
        ans <- ans & !pairedAns
      } else if(assertive::is_false(unpaired)){
        ans <- ans & pairedAns
      }
    } else {
      if(mismatches == TRUE){
        addAns <- vapply(strList,
                         function(str){
                           .is_continous_evenly_spaced(str[str$forward != 0]$forward) !=
                             .is_continous_evenly_spaced(str[str$reverse != 0]$reverse)
                         },
                         logical(1))
      }
      if(mismatches == FALSE){
        addAns <- vapply(strList,
                         function(str){
                           .is_continous_evenly_spaced(str[str$forward != 0]$forward) !=
                             .is_continous_evenly_spaced(str[str$reverse != 0]$reverse)
                         },
                         logical(1))
      }
      if(bulged == TRUE){
        addAns <- vapply(strList,
                         function(str){
                           length(str[str$forward != 0]$forward) !=
                             length(str[str$reverse != 0]$reverse)
                         },
                         logical(1))
      }
      if(bulged == FALSE){
        addAns <- vapply(strList,
                         function(str){
                           .is_continous_evenly_spaced(str[str$forward != 0]$forward) !=
                             .is_continous_evenly_spaced(str[str$reverse != 0]$reverse)
                         },
                         logical(1))
      }
      ans <- ans & addAns
    }
  }
  # 
  if(!is.na(length)){
    assertive::assert_all_are_whole_numbers(length)
    isLength <- lapply(strList,
                       function(str){
                         max(max(str[str$forward != 0,]$forward) - 
                               min(str[str$forward != 0,]$forward),
                             max(str[str$reverse != 0,]$reverse) - 
                               min(str[str$reverse != 0,]$reverse)) + 1
                       })
    ans <- ans & isLength == length
  }
  return(unname(ans))
}

# detect bulges in stem region and pad on opposite site
.has_unpaired_in_stem_region <- function(seqs,
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
