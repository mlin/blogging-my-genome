library(VariantAnnotation)
library(compiler)
library(data.table)
library(ggplot2)

# Load ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt
indivs <- read.delim("/20130606_sample_info.txt")
indivs <- indivs[,c('Sample','Population')]
indivs <- data.table(indivs)
setkey(indivs,Sample)

# Load 1000 Genomes VCF
keep <- list()
keep$fixed <- c("ALT")
keep$info <- c("AF","AA")
keep$geno <- c("GT")
kg <- readVcf("/pop.vcf.gz","hs37d5",ScanVcfParam(fixed=keep$fixed,info=keep$info,geno=keep$geno))
dim(kg)

# Filter 1000 Genomes VCF to biallelic SNVs with an inferred ancestral allele
is.nuc <- function(c) {
  C <- toupper(as.character(c))
  C == "A" | C == "G" | C == "C" | C == "T"
}
is.nuc <- cmpfun(is.nuc)
kg <- kg[elementLengths(info(kg)$AF) == 1,]
kg <- kg[sapply(ref(kg),is.nuc) & sapply(alt(kg),function(lst) { is.nuc(lst[[1]]) }),]
AA.known <- is.nuc(info(kg)$AA)
AA.known[is.na(AA.known)] <- FALSE
AA.upper <- toupper(info(kg)$AA)
kg <- kg[AA.known & (as.character(ref(kg)) == AA.upper | as.character(unlist(alt(kg))) == AA.upper),]
AA.upper <- toupper(info(kg)$AA)
dim(kg)

# Compute DAF
alt.is.ancestral <- as.character(unlist(alt(kg))) == AA.upper
DAF <- as.vector(info(kg)$AF)
DAF[alt.is.ancestral] <- 1-DAF[alt.is.ancestral]

# Select 5% < DAF < 50%
DAF.medium <- 0.05 <= DAF & DAF <= 0.5
kg <- kg[DAF.medium,]
alt.is.ancestral <- alt.is.ancestral[DAF.medium]
DAF <- DAF[DAF.medium]
AA.upper <- AA.upper[DAF.medium]
dim(kg)

# assign string positions such as "17:1234567" to facilitate intersection with focal VCF
vcf.positions <- function(vcf) {
  apply(cbind(as.character(seqnames(rowData(vcf))),start(ranges(rowData(vcf)))),1,function(x) { paste(x,collapse=":") } )
}
vcf.positions <- cmpfun(vcf.positions)
kg.positions <- vcf.positions(kg)

length(alt.is.ancestral)
length(kg.positions)
stopifnot(length(alt.is.ancestral) == nrow(kg))
stopifnot(length(kg.positions) == nrow(kg))

# Cons genotype matrix G with G_ij = # copies of derived allele j in individual i (0, 1, or 2)
make.G <- function(kg,alt.is.ancestral,positions) {
  GT <- geno(kg)$GT
  G <- matrix(1,nrow=ncol(GT),ncol=nrow(GT))
  rownames(G) <- colnames(GT)
  colnames(G) <- positions
  for (i in 1:nrow(G)) {
    G[i, GT[,i] == "0|0" & alt.is.ancestral] <- 2
    G[i, GT[,i] == "0|0" & !alt.is.ancestral] <- 0
    G[i, GT[,i] == "1|1" & alt.is.ancestral] <- 0
    G[i, GT[,i] == "1|1" & !alt.is.ancestral] <- 2
  }
  G
}
make.G <- cmpfun(make.G)
G <- make.G(kg,alt.is.ancestral,kg.positions)
dim(G)
# Vector of three-letter population identifiers corresponding to rows of G
pop <- indivs[rownames(G),][,Population]

# make a g vector (row of G) representing myself
me <- readVcf("/your.vcf.gz","hs37d5",ScanVcfParam(fixed=keep$fixed,info=keep$info,geno=keep$geno))
me.positions <- vcf.positions(me)
kg.in.me <- match(kg.positions,me.positions)
sum(!is.na(kg.in.me))
GT.in.me <- geno(me)$GT[kg.in.me]
# By default (for those positions with no entry in my VCF), fallaciously assume
# I'm homozygous for the reference allele.
me.g <- ifelse(alt.is.ancestral,2,0)
# Fill me.g where I do have calls.
me.g[!is.na(kg.in.me)] <- 1
me.g[!is.na(kg.in.me) & (GT.in.me == "0/0" | GT.in.me == "0|0") & alt.is.ancestral] <- 2
me.g[!is.na(kg.in.me) & (GT.in.me == "0/0" | GT.in.me == "0|0") & !alt.is.ancestral] <- 0
me.g[!is.na(kg.in.me) & (GT.in.me == "1/1" | GT.in.me == "1/1") & alt.is.ancestral] <- 0
me.g[!is.na(kg.in.me) & (GT.in.me == "1/1" | GT.in.me == "1/1") & !alt.is.ancestral] <- 2

# PCA on G
#pca.selector <- rep(TRUE,nrow(G))
pca.selector <- pop == 'CHB' | pop == 'CHS' | pop == 'JPT'
sum(pca.selector)
pca.varying <- sapply(1:ncol(G), function(j) { var(G[pca.selector,j]) > 0 })
sum(pca.varying)
pca <- prcomp(G[pca.selector,pca.varying],scale=FALSE,tol=0.1)

# place me on the PC axes
me.pca <- predict(pca,matrix(me.g[pca.varying],nrow=1,dimnames=list('YOUR_NAME',kg.positions[pca.varying])))

# Make some plots
pca.pop <- pop[pca.selector]
pca.pops <- sort(unique(pca.pop))
pd.shape <- numeric(length(pca.pop))
for(i in 1:length(pca.pops)) { pd.shape[pca.pop == pca.pops[i]] <- i-1 }
pd <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], PC3=pca$x[,3], population=pca.pop, shape=pd.shape)
popnames <- list(
  ASW='African-Amer',
  CEU='European',
  CHB='Chinese, Beijing',
  CHS='Chinese, Southern',
  CLM='Colombian',
  FIN='European',
  GBR='European',
  IBS='European',
  JPT='Japanese',
  MXL='Mexican-Amer',
  PUR='Puerto Ricans',
  TSI='European',
  YRI='African',
  LWK='African'
)
for (i in 1:length(popnames)) {
  levels(pd$population)[levels(pd$population) == names(popnames)[i]] <- popnames[[names(popnames)[i]]]
}
pd$population <- factor(pd$population,levels(pd$population)[c(13,2,5,7,18,15,6)])

png("/figures/unlabeled.png",width=700,height=512)
ggplot(pd, aes(PC1,PC2)) + geom_point() + theme(panel.background=element_blank())
dev.off()

cbbPalette <- c("#E69F00", "#CC79A7", "#D55E00", "#56B4E9", "#009E73", "#999999", "#0072B2", "#FF0000")
png("/figures/labeled.png",width=768,height=512)
ggplot(pd, aes(PC1,PC2,shape=population,color=population)) + geom_point() + scale_shape_manual(values=(1:nrow(pd))+1) + scale_color_manual(values=cbbPalette) + theme(panel.background=element_blank(), legend.text=element_text(size=10))
dev.off()

png("/figures/labeled_with_you.png",width=768,height=512)
ggplot(rbind(pd,data.frame(PC1=me.pca[,1],PC2=me.pca[,2],PC3=me.pca[,3],population='Mike Lin',shape=max(pd$shape)+1)),
       aes(PC1,PC2,shape=population,color=population)) + geom_point() + scale_shape_manual(values=(1:nrow(pd))+1) + scale_color_manual(values=cbbPalette) + theme(panel.background=element_blank(), legend.text=element_text(size=10))
dev.off()
