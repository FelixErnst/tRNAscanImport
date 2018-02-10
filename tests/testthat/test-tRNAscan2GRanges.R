library(tRNAscanImport)

test_that("tests:",{
  file <- system.file("extdata", 
                      file = "sacCer3-tRNAs.ss.sort", 
                      package = "tRNAscanImport")
  lines <- readLines(con = file, n = 7L)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines)
  expected <- c("trna","type","intron", "pseudogene", "hmm", "secstruct", 
                "infernal", "seq","str")
  expect_named(actual, expected)
  expect_equal(min(unlist(lapply(actual, length))), 0)
  expect_equal(max(unlist(lapply(actual, length))), 6)
  
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines[c(1:3,5:7)])
  expected <- c("trna","type","intron","seq","str")
  expect_named(actual, expected)
  
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines[c(1:2,4:7)])
  expected <- c("trna","type","pseudogene", "hmm", "secstruct", 
                "infernal","seq","str")
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
  
  df <- tRNAscanImport:::.read_tRNAscan(file)
  actual <- tRNAscanImport:::.cut_introns(df)
  expect_false(identical(df[1,"tRNA_seq"], actual[1,"tRNA_seq"]))
  expect_true(identical("GGGCGTGTGGTCTAGTGGTATGATTCTCGCTTTGGGTGCGAGAGGcCCTGGGTTCAATTCCCAGCTCGCCCC", 
                        actual[1,"tRNA_seq"]))
  expect_false(identical(df[1,"seq"], actual[1,"tRNA_str"]))
  expect_true(identical(">>>>>.>..>>>.........<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<.<<<<<.", 
                        actual[1,"tRNA_str"]))
  
  gr <- tRNAscanImport:::tRNAscanImport(file)
  length <- as.numeric(S4Vectors::mcols(gr)$tRNA_length)
  intron_locstart <- as.numeric(S4Vectors::mcols(gr)$tRNAscan_intron.locstart)
  intron_locstart[is.na(intron_locstart)] <- 0
  intron_locend <- as.numeric(S4Vectors::mcols(gr)$tRNAscan_intron.locend)
  intron_locend[is.na(intron_locend)] <- 0
  intron_length <- intron_locend - intron_locstart
  length <- length - vapply(intron_length, function(l){
    if(l > 0){
      return(l + 1)
    }
    0
  },numeric(1))
  expect_equal(length,BiocGenerics::width(S4Vectors::mcols(gr)$tRNA_seq))
  expect_equal(length,BiocGenerics::width(S4Vectors::mcols(gr)$tRNA_str))
})

test_that("input failure test:",{
  expect_error(
    tRNAscanImport:::.parse_tRNAscan_block(),
    'argument "lines" is missing'
  )
  
  file <- tempfile()
  writeLines(c(""),file)
  expect_error(
    actual <- tRNAscanImport:::.parse_tRNAscan(file),
    'Empty file.'
  )
  file <- tempfile()
  writeLines(c("\n\n\n\n"),file)
  expect_error(
    actual <- tRNAscanImport:::.parse_tRNAscan(file),
    'No tRNA information detected. Please make sure'
  )
  
  file <- system.file("extdata", 
                      file = "sacCer3-tRNAs.ss.sort", 
                      package = "tRNAscanImport")
  lines <- readLines(con = file, n = 7L)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines[1:2])
  expect_named(actual, NULL)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines[c(1:2,5:7)])
  expect_named(actual, NULL)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(list())
  expect_named(actual, NULL)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(c())
  expect_named(actual, NULL)
  
  expect_error(
    tRNAscanImport:::.has_CCA_end(),
    'argument "seq" is missing'
  )
  expect_false(
    tRNAscanImport:::.has_CCA_end("")
  )
  expect_error(
    tRNAscanImport:::.read_tRNAscan(),
    'argument "file" is missing'
  )
  expect_error(
    tRNAscanImport:::.cut_introns(),
    'argument "df" is missing'
  )
})