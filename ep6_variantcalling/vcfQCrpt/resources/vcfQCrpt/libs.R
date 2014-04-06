local({r <- getOption("repos");  r["CRAN"] <- "http://cran.rstudio.org"; options(repos=r)})
install.packages(c("knitr","markdown","ggplot2","RCurl","bitops"),quiet=TRUE)
source("http://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
