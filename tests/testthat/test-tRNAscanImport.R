library(tRNAscanImport)

context("tests")
test_that("tests:",{
  file <- system.file("extdata", 
                      file = "yeast.tRNAscan", 
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
  expect_true(identical(
    "GGGCGTGTGGTCTAGTGGTATGATTCTCGCTTTGGGTGCGAGAGGcCCTGGGTTCAATTCCCAGCTCGCCCC", 
    actual[1,"tRNA_seq"]))
  expect_false(identical(df[1,"seq"], actual[1,"tRNA_str"]))
  expect_true(identical(
    ">>>>>.>..>>>.........<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<.<<<<<.", 
    actual[1,"tRNA_str"]))
  
  gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  length <- as.numeric(S4Vectors::mcols(gr)$tRNA_length)
  expect_equal(length,BiocGenerics::width(S4Vectors::mcols(gr)$tRNA_seq))
  expect_equal(length,BiocGenerics::width(S4Vectors::mcols(gr)$tRNA_str))
})

context("type tests gr")
test_that("type tests gr:",{
  file <- system.file("extdata", 
                      file = "yeast.tRNAscan", 
                      package = "tRNAscanImport")
  gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  expect_true(istRNAscanGRanges(gr))
  expect_type(mcols(gr)$no, "integer")
  expect_type(mcols(gr)$tRNA_length, "integer")
  expect_type(mcols(gr)$tRNA_type, "character")
  expect_type(mcols(gr)$tRNA_anticodon, "character")
  expect_type(mcols(gr)$tRNA_anticodon.start, "integer")
  expect_type(mcols(gr)$tRNA_anticodon.end, "integer")
  expect_type(mcols(gr)$tRNAscan_score, "double")
  expect_type(mcols(gr)$tRNA_seq, "S4")
  expect_type(mcols(gr)$tRNA_str, "S4")
  expect_type(mcols(gr)$tRNA_CCA.end, "logical")
  expect_type(mcols(gr)$tRNAscan_potential.pseudogene, "logical")
  expect_type(mcols(gr)$tRNAscan_intron.start, "integer")
  expect_type(mcols(gr)$tRNAscan_intron.end, "integer")
  expect_type(mcols(gr)$tRNAscan_intron.locstart, "integer")
  expect_type(mcols(gr)$tRNAscan_intron.locend, "integer")
  expect_type(mcols(gr)$tRNAscan_hmm.score, "double")
  expect_type(mcols(gr)$tRNAscan_sec.str.score, "double")
  expect_type(mcols(gr)$tRNAscan_infernal, "double")
  expect_true(all(gr$tRNAscan_potential.pseudogene == FALSE))
})

context("type tests gff")
test_that("type tests gff:",{
  file <- system.file("extdata", 
                      file = "yeast.tRNAscan", 
                      package = "tRNAscanImport")
  gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  gff <- tRNAscanImport::tRNAscan2GFF(gr)
  expect_type(mcols(gff)$source, "integer")
  expect_true(is.factor(mcols(gff)$source))
  expect_type(mcols(gff)$type, "integer")
  expect_true(is.factor(mcols(gff)$type))
  expect_type(mcols(gff)$score, "integer")
  expect_type(mcols(gff)$phase, "integer")
  expect_type(mcols(gff)$ID, "character")
  # 
  expect_type(mcols(gff)$no, "integer")
  expect_type(mcols(gff)$tRNA_length, "integer")
  expect_type(mcols(gff)$tRNA_type, "character")
  expect_type(mcols(gff)$tRNA_anticodon, "character")
  expect_type(mcols(gff)$tRNA_anticodon.start, "integer")
  expect_type(mcols(gff)$tRNA_anticodon.end, "integer")
  expect_type(mcols(gff)$tRNAscan_score, "double")
  expect_type(mcols(gff)$tRNA_seq, "character")
  expect_type(mcols(gff)$tRNA_str, "character")
  expect_type(mcols(gff)$tRNA_CCA.end, "logical")
  expect_type(mcols(gff)$tRNAscan_potential.pseudogene, "logical")
  expect_type(mcols(gff)$tRNAscan_intron.start, "integer")
  expect_type(mcols(gff)$tRNAscan_intron.end, "integer")
  expect_type(mcols(gff)$tRNAscan_intron.locstart, "integer")
  expect_type(mcols(gff)$tRNAscan_intron.locend, "integer")
  expect_type(mcols(gff)$tRNAscan_hmm.score, "double")
  expect_type(mcols(gff)$tRNAscan_sec.str.score, "double")
  expect_type(mcols(gff)$tRNAscan_infernal, "double")
})

context("input failure tests")
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
                      file = "yeast.tRNAscan", 
                      package = "tRNAscanImport")
  lines <- readLines(con = file, n = 7L)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines[1:2])
  expect_named(actual, NULL)
  actual <- tRNAscanImport:::.parse_tRNAscan_block(lines[c(1:2,6:7)])
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
