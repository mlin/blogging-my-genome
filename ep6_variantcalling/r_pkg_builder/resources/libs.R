local({r <- getOption("repos");  r["CRAN"] <- "http://cran.rstudio.org"; options(repos=r)})
cran.pkgs <- as.character(read.delim("/cran.txt",header=FALSE)[,1])
if (length(cran.pkgs) > 0) {
	install.packages(cran.pkgs)
}
bioc.pkgs <- as.character(read.delim("/bioconductor.txt",header=FALSE)[,1])
if (length(bioc.pkgs) > 0) {
	source("http://bioconductor.org/biocLite.R")
	biocLite(bioc.pkgs)
}
