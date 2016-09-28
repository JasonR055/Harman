# Harman

The removal of batch effects from datasets using a PCA and constrained
optimisation based technique.

Harman is a PCA and constrained optimisation based technique that maximises the
removal of batch effects from datasets, with the constraint that the probability
of overcorrection (i.e. removing genuine biological signal along with batch
noise) is kept to a fraction which is set by the end-user (Oytam et al, 2016; DOI:10.1186/s12859-016-1212-5).

The Harman approach can be generalised to any high dimensional dataset where an
experimental factor of interest to be kept is declared and another (usually
technical) factor is declared to be removed, with the constraint that the factor
to be kept and removed are not completely confounded.

To compile the package on a command line:  
__Windows__: `Rcmd.exe INSTALL --preclean --no-multiarch --with-keep.source Harman`  
__Linux__: `R CMD INSTALL --preclean --build Harman -l ~/R/library`  

To set up help files:  
`roxygen2::roxygenize(package.dir = "/path/to/harman")`  

_N.B._  
* The `roxygenize()` function should create the `NAMESPACE` file also.  

To use harman:
```R
install.packages("Harman")
library("Harman")
data(OLF)
olf.harman <- harman(olf.data, expt=olf.info$Treatment, batch=olf.info$Batch)
olf.data.corrected <- reconstructData(olf.harman)
```
 
Required packages for `R CMD Check` to pass (Packages used in Vignettes)
```R
source("http://www.bioconductor.org/biocLite.R")
biocLite()
biocLite(c("BiocGenerics", "BiocStyle", "knitr", "rmarkdown", "RUnit",
"missMethyl", "RColorBrewer", "bladderbatch", "limma", "minfi", "lumi",
"msmsEDA", "affydata", "minfiData"))
```
Also required to build vignettes is a version of Texlive for the pdflatex
program and a version of the NetCDF library.

**Further documentation and the Matlab version is available here:** http://www.bioinformatics.csiro.au/harman/