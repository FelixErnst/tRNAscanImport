#' @include tRNA.R
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
#' @param length the length as integer
#' @param unpaired logical: has unpaired nucleotides
#' @param paired logical: has paired nucleotides (only used for loops)
#' @param mismatches logical: has mismatched nucleotides 
#' @param bulged logical: has mismatched nucleotides of different length
#' creating a bulge
#' 
#' @return a logical vector of the length or input GRanges object
#' 
#' @examples
#' gr <- import.tRNAscanAsGRanges(system.file("extdata", 
#'                                file = "sacCer3-tRNAs.ss.sort", 
#'                                package = "tRNAscanImport"))
#' hasTStem(gr, length = 5, mismatches = TRUE)
#' gr[hasTStem(gr, length = 5, mismatches = TRUE)]
NULL

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasTStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("TStem",
                                                  gr,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasDStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("DStem",
                                                  gr,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasAcceptorStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("acceptorStem",
                                                  gr,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

#' @rdname tRNA-subset
#' 
#' @export
setMethod(
  f = "hasAnticodonStem",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length,
                        unpaired,
                        mismatches,
                        bulged) .subset_tRNA_stem("anticodonStem",
                                                  gr,
                                                  length,
                                                  unpaired,
                                                  mismatches,
                                                  bulged))

.subset_tRNA_stem <- function(ident,
                              gr,
                              length,
                              unpaired,
                              mismatches,
                              bulged){
  ans <- rep(TRUE,length(gr))
  # input check
  .check_trnascan_granges(gr, TRNA_FEATURES)
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
  if(!is.na(unpaired) && 
     is.na(mismatches) && 
     is.na(bulged)){
    assertive::assert_is_a_bool(unpaired)
    if(assertive::is_true(unpaired)){
      mismatches <- TRUE
      bulged <- TRUE
    } else if(assertive::is_false(unpaired)) {
      mismatches <- FALSE
      bulged <- FALSE
    }
  }
  # get structure information
  strList <- .get_base_pairing(gr$tRNA_str)
  str <- .get_tRNA_structures(tRNAStructureFunctionList[[ident]],
                              gr,
                              strList)[[1]]
  if(!is.na(unpaired) ||
     !is.na(length)){
    strList <- mapply(.subset_structure,
                      split(c(str$prime5,str$prime3),
                            c(1:length(str$prime5),1:length(str$prime3))),
                      strList,
                      SIMPLIFY = FALSE)
  }
  # apply structure subsetting
  if(!is.na(unpaired)){
    #
    forwardContinously <- vapply(
      strList,
      function(str){
        .is_continous_evenly_spaced(str[str$forward != 0,]$forward)
      },
      logical(1))
    reverseContinously <- vapply(
      strList,
      function(str){
        .is_continous_evenly_spaced(str[str$reverse != 0,]$reverse)
      },
      logical(1))
    forwardLength <- vapply(
      strList,
      function(str){
        max(str[str$forward != 0,]$forward) - 
          min(str[str$forward != 0,]$forward) + 1
      },
      double(1))
    reverseLength <- vapply(
      strList,
      function(str){
        max(str[str$reverse != 0,]$reverse) - 
          min(str[str$reverse != 0,]$reverse) + 1
      },
      double(1))
    #
    ansMismatches <- ans
    ansBulged <- ans
    if(assertive::is_false(mismatches)){
      ansMismatches <- forwardLength == reverseLength &
        forwardContinously &
        reverseContinously
    } else if(assertive::is_true(mismatches)) {
      ansMismatches <- forwardLength == reverseLength &
        ((forwardContinously & !reverseContinously) | 
           (!forwardContinously & reverseContinously) |
           (!forwardContinously & !reverseContinously))
    }
    if(assertive::is_false(bulged)){
      ansBulged <- forwardLength == reverseLength &
        ((forwardContinously & reverseContinously) | 
           (!forwardContinously & !reverseContinously))
    } else if(assertive::is_true(bulged)) {
      ansBulged <- forwardLength != reverseLength &
        ((forwardContinously & !reverseContinously) | 
           (!forwardContinously & reverseContinously) |
           (!forwardContinously & !reverseContinously))
    }
    if(assertive::is_true(mismatches) && 
       assertive::is_true(bulged)){
      ans <- ans & (ansMismatches | ansBulged)
    } else {
      ans <- ans & ansMismatches & ansBulged
    }
  }
  # apply length subsetting 
  if(!is.na(length)){
    assertive::assert_all_are_whole_numbers(length)
    isLength <- lapply(strList,
                       function(str){
                         max(max(str[str$forward != 0,]$forward) - 
                               min(str[str$forward != 0,]$forward),
                             max(str[str$reverse != 0,]$reverse) - 
                               min(str[str$reverse != 0,]$reverse)) + 1
                       })
    ans <- ans & 
      isLength == length
  }
  return(unname(ans))
}


# loop subsetting --------------------------------------------------------------

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasTloop",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length) .subset_tRNA_loop("Tloop",
                                                  gr,
                                                  length,
                                                  NA,
                                                  NA,
                                                  NA))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasDloop",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length) .subset_tRNA_loop("Dloop",
                                                  gr,
                                                  length,
                                                  NA,
                                                  NA,
                                                  NA))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasAnticodonLoop",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length) .subset_tRNA_loop("anticodonLoop",
                                                  gr,
                                                  length,
                                                  NA,
                                                  NA,
                                                  NA))

#' @rdname tRNA-subset
#' @export
setMethod(
  f = "hasVariableLoop",
  signature = signature(gr = "GRanges"),
  definition = function(gr,
                        length,
                        paired,
                        mismatches,
                        bulged) .subset_tRNA_loop("variableLoop",
                                                  gr,
                                                  length,
                                                  paired,
                                                  mismatches,
                                                  bulged))

.subset_tRNA_loop <- function(ident,
                              gr,
                              length,
                              paired,
                              mismatches,
                              bulged){
  ans <- rep(TRUE,length(gr))
  # input check
  .check_trnascan_granges(gr, TRNA_FEATURES)
  if(!(ident %in% names(tRNAStructureFunctionList))){
    stop("Unknown identifier '",ident,"'.")
  }
  if(!is.na(mismatches)){
    assertive::assert_is_a_bool(mismatches)
    paired <- TRUE
  }
  if(!is.na(bulged)){
    assertive::assert_is_a_bool(bulged)
    paired <- TRUE
  }
  if(!is.na(paired)){
    assertive::assert_is_a_bool(paired)
  }
  #
  strList <- .get_base_pairing(gr$tRNA_str)
  str <- .get_tRNA_structures(tRNAStructureFunctionList[[ident]],
                              gr,
                              strList)[[1]]
  if(!is.na(paired) ||
     !is.na(length)){
    strList <- mapply(.subset_structure,
                      split(str,1:length(str)),
                      strList,
                      MoreArgs = list(pairedOnly = FALSE),
                      SIMPLIFY = FALSE)
  }
  # apply paired subsetting
  if(!is.na(paired)){
    pairedAns <- vapply(
      strList,
      function(str){
        any(str$forward > 0)
      },
      logical(1))
    # if any paired nucleotides are found
    if(any(pairedAns) &&
       (!is.na(mismatches) ||
        !is.na(bulged))){
      forwardContinously <- vapply(
        strList,
        function(str){
          .is_continous_evenly_spaced(str[str$forward != 0,]$forward)
        },
        logical(1))
      reverseContinously <- vapply(
        strList,
        function(str){
          .is_continous_evenly_spaced(str[str$reverse != 0,]$reverse)
        },
        logical(1))
      forwardLength <- vapply(
        strList,
        function(str){
          max(str[str$forward != 0,]$forward) - 
            min(str[str$forward != 0,]$forward) + 1
        },
        double(1))
      reverseLength <- vapply(
        strList,
        function(str){
          max(str[str$reverse != 0,]$reverse) - 
            min(str[str$reverse != 0,]$reverse) + 1
        },
        double(1))
      #
      ansMismatches <- ans
      ansBulged <- ans
      if(assertive::is_false(mismatches)){
        ansMismatches <- forwardLength == reverseLength &
          forwardContinously &
          reverseContinously
      } else if(assertive::is_true(mismatches)) {
        ansMismatches <- forwardLength == reverseLength &
          ((forwardContinously & !reverseContinously) | 
             (!forwardContinously & reverseContinously) |
             (!forwardContinously & !reverseContinously))
      }
      if(assertive::is_false(bulged)){
        ansBulged <- forwardLength == reverseLength &
          ((forwardContinously & reverseContinously) | 
             (!forwardContinously & !reverseContinously))
      } else if(assertive::is_true(bulged)) {
        ansBulged <- forwardLength != reverseLength &
          ((forwardContinously & !reverseContinously) | 
             (!forwardContinously & reverseContinously) |
             (!forwardContinously & !reverseContinously))
      }
      if(assertive::is_true(mismatches) && 
         assertive::is_true(bulged)){
        ans <- ans & (ansMismatches | ansBulged)
      } else {
        ans <- ans & ansMismatches & ansBulged
      }
    }
    ans <- ans & 
      pairedAns
  }
  # apply length subsetting 
  if(!is.na(length)){
    assertive::assert_all_are_whole_numbers(length)
    isLength <- lapply(strList,
                       function(str){
                         max(str$pos) - 
                           min(str$pos) + 1
                       })
    ans <- ans & 
      isLength == length
  }
  return(unname(ans))
}