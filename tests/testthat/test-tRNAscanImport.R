library(tRNAscanImport)

context("tests")
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
  
  gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
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

context("type tests gr")
test_that("type tests gr:",{
  file <- system.file("extdata", 
                      file = "sacCer3-tRNAs.ss.sort", 
                      package = "tRNAscanImport")
  gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  expect_type(mcols(gr)$no, "integer")
  expect_type(mcols(gr)$tRNA_length, "integer")
  expect_type(mcols(gr)$tRNA_type, "character")
  expect_type(mcols(gr)$tRNA_anticodon, "character")
  expect_type(mcols(gr)$tRNA_anticodon.start, "integer")
  expect_type(mcols(gr)$tRNA_anticodon.end, "integer")
  expect_type(mcols(gr)$tRNAscan_score, "double")
  expect_type(mcols(gr)$tRNA_seq, "S4")
  expect_type(mcols(gr)$tRNA_str, "character")
  expect_type(mcols(gr)$tRNA_CCA.end, "logical")
  expect_type(mcols(gr)$tRNAscan_potential.pseudogene, "logical")
  expect_type(mcols(gr)$tRNAscan_intron.start, "integer")
  expect_type(mcols(gr)$tRNAscan_intron.end, "integer")
  expect_type(mcols(gr)$tRNAscan_intron.locstart, "integer")
  expect_type(mcols(gr)$tRNAscan_intron.locend, "integer")
  expect_type(mcols(gr)$tRNAscan_hmm.score, "double")
  expect_type(mcols(gr)$tRNAscan_sec.str.score, "double")
  expect_type(mcols(gr)$tRNAscan_infernal, "double")
})

context("type tests gff")
test_that("type tests gff:",{
  file <- system.file("extdata", 
                      file = "sacCer3-tRNAs.ss.sort", 
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
                      file = "sacCer3-tRNAs.ss.sort", 
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

context("tRNA structure seqs")
test_that("tRNA structure seqs:",{
  file <- system.file("extdata", 
                      file = "sacCer3-tRNAs.ss.sort", 
                      package = "tRNAscanImport")
  tRNA <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  tRNA <- tRNA[1]
  actual <- tRNAscanImport:::.getAcceptorStem(tRNA)
  expect_named(actual, c("prime5","prime3"))
  expect_named(actual[[1]], c("start","end"))
  expect_named(actual[[2]], c("start","end"))
  expect_equal(actual[[2]], list(start = 65,end = 71))
  actual <- tRNAscanImport:::.getDiscriminator(tRNA)
  expect_equal(actual, 72)
  actual <- tRNAscanImport:::.getDstem(tRNA)
  expect_named(actual, c("prime5","prime3"))
  expect_named(actual[[1]], c("start","end"))
  expect_named(actual[[2]], c("start","end"))
  expect_equal(actual[[2]], list(start = 22,end = 25))
  actual <- tRNAscanImport:::.getAnticodonStem(tRNA)
  expect_named(actual, c("prime5","prime3"))
  expect_named(actual[[1]], c("start","end"))
  expect_named(actual[[2]], c("start","end"))
  expect_equal(actual[[2]], list(start = 38,end = 42))
  actual <- tRNAscanImport:::.getTstem(tRNA)
  expect_named(actual, c("prime5","prime3"))
  expect_named(actual[[1]], c("start","end"))
  expect_named(actual[[2]], c("start","end"))
  expect_equal(actual[[2]], list(start = 60,end = 64))
  actual <- tRNAscanImport:::.getDloop(tRNA)
  expect_named(actual, c("start","end"))
  expect_equal(actual, list(start = 13,end = 21))
  actual <- tRNAscanImport:::.getAnticodonloop(tRNA)
  expect_named(actual, c("start","end"))
  expect_equal(actual, list(start = 31,end = 37))
  actual <- tRNAscanImport:::.getVariableLoop(tRNA)
  expect_named(actual, c("start","end"))
  expect_equal(actual, list(start = 43,end = 47))
  actual <- tRNAscanImport:::.getTloop(tRNA)
  expect_named(actual, c("start","end"))
  expect_equal(actual, list(start = 53,end = 59))
  
  tRNA <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  seqs <- tRNAscanImport::gettRNAstructureSeqs(tRNA,
                                               joinCompletely = TRUE,
                                               joinFeatures = FALSE)
  seqsChars <- gsub("-","",as.character(seqs))
  expect_equal(seqsChars,as.character(tRNA$tRNA_seq))
})

