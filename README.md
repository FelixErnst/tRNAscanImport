# tRNAscanImport

The default tRNAscan-SE ([Lowe et al. 1997](#Literature)) output is formatted
text document containing text blocks per tRNA delimited by an empty line. To
access the information in a BioC context the conversion to a GRanges object
comes to mind. This task is performed by `import.tRNAscanAsGRanges()`, which
uses regular expressions to extract the information from the text blocks. The
result can be used directly or saved as gff3 file for further use.
Please refer to the vignette for an example usage case.

# Installation

The release version can be installed using `biocLite` from bioconductor.

```r
source("https://bioconductor.org/biocLite.R")
biocLite("tRNAscanImport")
```

The development version is of course accessible via GitHub and the `devtools`
package
```r
devtools::install_github("FelixErnst/tRNAscanImport")
```

# Literature

Depending on the development on tRNAscan-SE this might become redundant, since
a gff3 export by tRNAscan-SE might resolve the conversion issue. 

- Lowe, T.M. & Eddy, S.R. (1997) tRNAscan-SE: A program for 
  improved detection of transfer RNA genes in genomic sequence. Nucl. Acids Res. 
  25: 955-964.