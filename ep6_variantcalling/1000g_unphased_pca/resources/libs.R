local({r <- getOption("repos");  r["CRAN"] <- "http://cran.rstudio.org"; options(repos=r)})
# bitops & RCurl needed for VariantAnnotation
install.packages(c("ggplot2", "data.table", "bitops", "RCurl"),quiet=TRUE)
source("http://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
