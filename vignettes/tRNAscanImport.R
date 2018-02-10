## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown(css.files = c('custom.css'))

## ------------------------------------------------------------------------
library(tRNAscanImport)
sacCer3_file <- system.file("extdata", 
                            file = "sacCer3-tRNAs.ss.sort", 
                            package = "tRNAscan2GRanges")

# output for sacCer3
# Before
readLines(con = sacCer3_file, n = 7L)

## ------------------------------------------------------------------------
# output for sacCer3
# After
gr <- import.tRNAscanAsGRanges(sacCer3_file)
head(gr, 2)

## ------------------------------------------------------------------------
library(Biostrings, quietly = TRUE)
library(rtracklayer, quietly = TRUE)
# suppressMessages(library(rtracklayer, quietly = TRUE))
# Save tRNA sequences
writeXStringSet(gr$tRNA_seq, filepath = tempfile())
# to be GFF3 compliant use tRNAscan2GFF
gff <- tRNAscan2GFF(gr)
export.gff3(gff, con = tempfile())

## ------------------------------------------------------------------------
library(GenomicRanges, quietly = TRUE)
# tRNAscan-SE output for hg38
hg38_file <- system.file("extdata", 
                         file = "hg38-tRNAs.ss.sort", 
                         package = "tRNAscan2GRanges")
# tRNAscan-SE output for E. coli MG1655
eco_file <- system.file("extdata", 
                        file = "eschColi_K_12_MG1655-tRNAs.ss.sort", 
                        package = "tRNAscan2GRanges")
# import tRNAscan-SE files
gr_hg <- import.tRNAscanAsGRanges(hg38_file)
gr_eco <- import.tRNAscanAsGRanges(eco_file)

# get summary plots if ggplot2 is installed
grl <- GRangesList(Sce = gr, 
                   Hsa = gr_hg, 
                   Eco = gr_eco)
plots <- gettRNAscanPlots(grl)

## ---- fig.cap = "tRNA length."-------------------------------------------
plots$length

## ---- fig.cap = "tRNAscan-SE scores."------------------------------------
plots$tRNAscan_score

## ---- fig.cap = "tRNA GC content."---------------------------------------
plots$gc

## ---- fig.cap = "tRNAs with introns."------------------------------------
plots$introns

## ------------------------------------------------------------------------
sessionInfo()

