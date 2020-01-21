library(tRNAscanImport)

context("tRNAscan import")
test_that("tRNAscan import:",{
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

context("tRNAscan data integrity")
test_that("tRNAscan data integrity:",{
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
  #
  mcols(gr)$tRNA_seq <- NULL
  expect_warning(expect_false(istRNAscanGRanges(gr)),
                 "Input GRanges object does not meet the requirements")
  expect_error(istRNAscanGRanges(data.frame()))
  expect_warning(expect_false(.check_trnascan_granges(data.frame(),"")),
                 "Input is not a GRanges object")
  #
  expect_warning(tRNAscanImport::import.tRNAscanAsGRanges(file, as.GFF3 = 1),
                 "'as.GFF3' is not a bool. Resetting")
  expect_warning(tRNAscanImport::import.tRNAscanAsGRanges(file,
                                                          trim.intron = 1),
                 "'trim.intron' is not a bool. Resetting")
})

context("tRNAscan gff integrity")
test_that("tRNAscan gff integrity:",{
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
  #
  gff2 <- tRNAscanImport::import.tRNAscanAsGRanges(file, as.GFF3 = TRUE)
  expect_equal(gff, gff2)
})

context("tRNA precursor")
test_that("tRNA precursor:",{
  library(BSgenome.Scerevisiae.UCSC.sacCer3)
  file <- system.file("extdata",
                      file = "yeast.tRNAscan",
                      package = "tRNAscanImport")
  gr <- tRNAscanImport::import.tRNAscanAsGRanges(file)
  genome <- getSeq(BSgenome.Scerevisiae.UCSC.sacCer3)
  #
  expect_error(get.tRNAprecursor(gr, genome),
               "Not all chromosomes referenced in the tRNAscan data are")
  #
  names(genome) <- c(names(genome)[-17],"chrmt")
  expect_equal(BSgenome.Scerevisiae.UCSC.sacCer3,
               tRNAscanImport:::.norm_genome(BSgenome.Scerevisiae.UCSC.sacCer3))
  expect_equal(genome,tRNAscanImport:::.norm_genome(genome))
  fl <- system.file("extdata", "ce2dict1.fa", package="Rsamtools",
                    mustWork=TRUE)
  ff <- Rsamtools::FaFile(fl)
  expect_equal(ff,tRNAscanImport:::.norm_genome(ff))
  #
  expect_error(get.tRNAprecursor(gr, genome, add.5prime = 1.2),
               "'add.5prime' and 'add.3prime' must be integer values")
  expect_error(get.tRNAprecursor(gr, genome, add.3prime = 1.2),
               "'add.5prime' and 'add.3prime' must be integer values")
  expect_error(get.tRNAprecursor(gr, genome, add.3prime = c(1,2)),
               "'add.5prime' and 'add.3prime' must be integer values")
  expect_error(get.tRNAprecursor(gr, genome, add.3prime = c(1L,2L)),
               "'add.3prime' must be a single integer value")
  expect_error(get.tRNAprecursor(gr, genome, add.5prime = c(1L,2L)),
               "'add.5prime' must be a single integer value")
  expect_warning(get.tRNAprecursor(gr, genome, trim.intron = 1),
                 "'trim.intron' is not a bool. Resetting ")
  #
  actual <- get.tRNAprecursor(gr, genome)
  genome2 <- BSgenome.Scerevisiae.UCSC.sacCer3
  seqnames(genome2) <- c(seqnames(genome2)[-17],"chrmt")
  actual2 <- get.tRNAprecursor(gr, genome2)
  expect_equal(as.character(actual), as.character(actual2))
  expect_equal(width(head(actual,n=5)), c(203L, 173L, 214L, 182L, 184L))
  #
  actual <- get.tRNAprecursor(gr, genome, trim.intron = TRUE)
  expect_equal(width(head(actual,n=5)), c(172L, 173L, 182L, 182L, 184L))
  actual <- get.tRNAprecursor(gr, genome, trim.intron = TRUE, add.5prime = 10L)
  expect_equal(width(head(actual,n=5)), c(92L, 93L, 102L, 102L, 104L))
  actual <- get.tRNAprecursor(gr, genome, trim.intron = TRUE, add.5prime = 10L,
                              add.3prime = 50L)
  expect_equal(width(head(actual,n=5)), c(132L, 133L, 142L, 142L, 144L))
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
  expect_error(tRNAscanImport:::.norm_genome(),
               'argument "genome" is missing, with no default')
  expect_error(tRNAscanImport:::.norm_genome(data.frame()),
               "'genome' must be an object of class 'BSgenome'")
})
