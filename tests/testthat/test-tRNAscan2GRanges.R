library(tRNAscan2GRanges)

test_that("tests:",{
  file <- system.file("extdata", 
                      file = "tRNAscan.sort", 
                      package = "tRNAscan2GRanges")
  lines <- readLines(con = file, n = 7L)
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(lines)
  expected <- c("trna","type","intron","hmm","seq","str")
  expect_named(actual, expected)
  expect_equal(min(unlist(lapply(actual, length))), 2)
  
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(lines[c(1:3,5:7)])
  expected <- c("trna","type","intron","seq","str")
  expect_named(actual, expected)
  
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(lines[c(1:2,4:7)])
  expected <- c("trna","type","hmm","seq","str")
  expect_named(actual, expected)
  
  expect_equal(.has_CCA_end(actual$seq[2],
                            actual$str[2]), 
               FALSE)
  expect_equal(.has_CCA_end(paste0(actual$seq[2],"CCA"),
                            paste0(actual$str[2])), 
               FALSE)
  expect_equal(.has_CCA_end(paste0(actual$seq[2]),
                            paste0(actual$str[2],"...")), 
               FALSE)
  expect_equal(.has_CCA_end(paste0(actual$seq[2],"CCA"),
                            paste0(actual$str[2],"...")), 
               TRUE)
  
  df <- tRNAscan2GRanges:::.read_tRNAscan(file)
  actual <- tRNAscan2GRanges:::.cut_introns(df)
  expect_false(identical(df[1,"seq"], actual[1,"seq"]))
  expect_true(identical("GGGCGTGTGGTCTAGTGGTATGATTCTCGCTTTGGGTGCGAGAGGcCCTGGGTTCAATTCCCAGCTCGCCCC", 
                        actual[1,"seq"]))
  expect_false(identical(df[1,"seq"], actual[1,"str"]))
  expect_true(identical(">>>>>.>..>>>.........<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<.<<<<<.", 
                        actual[1,"str"]))
})

test_that("input failure test:",{
  expect_error(
    tRNAscan2GRanges:::.parse_tRNAscan_block(),
    'argument "lines" is missing'
  )
  
  file <- tempfile()
  writeLines(c(""),file)
  expect_error(
    actual <- tRNAscan2GRanges:::.parse_tRNAscan(file),
    'Empty file.'
  )
  file <- tempfile()
  writeLines(c("\n\n\n\n"),file)
  expect_error(
    actual <- tRNAscan2GRanges:::.parse_tRNAscan(file),
    'No tRNA information detected. Please make sure'
  )
  
  file <- system.file("extdata", 
                      file = "tRNAscan.sort", 
                      package = "tRNAscan2GRanges")
  lines <- readLines(con = file, n = 7L)
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(lines[1:2])
  expect_named(actual, NULL)
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(lines[c(1:2,5:7)])
  expect_named(actual, NULL)
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(list())
  expect_named(actual, NULL)
  actual <- tRNAscan2GRanges:::.parse_tRNAscan_block(c())
  expect_named(actual, NULL)
  
  expect_error(
    tRNAscan2GRanges:::.has_CCA_end(),
    'argument "seq" is missing'
  )
  expect_false(
    tRNAscan2GRanges:::.has_CCA_end("")
  )
  expect_error(
    tRNAscan2GRanges:::.read_tRNAscan(),
    'argument "file" is missing'
  )
  expect_error(
    tRNAscan2GRanges:::.cut_introns(),
    'argument "df" is missing'
  )
})