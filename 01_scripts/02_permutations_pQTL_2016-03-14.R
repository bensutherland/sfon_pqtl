# Analysis of the phenotypic QTL from Sauvage et al with the new map from Sutherland et al.
# analysis v2.0
# Part 02/03 (permutations on cluster computer)
# in R/QTL (cp cross type)

# rm(list=ls())

# install libraries
require(qtl)
install.packages("rlecuyer")
require(rlecuyer) #random number generator in parallel
require(plyr)

# set working directory for the cluster computer
setwd("~/Documents/sfqtl_macpro_2016-09-16/sfon_pqtl")

# Load Part 1 results:
load("02_data/sfon_01_output.RData")

# Set permutation variables (perm should be 1000):
num.perm = 1000
num.cluster = 20

####### 2A Single QTL, no covariate #####
# Create scanone object for phenos not requiring covariate (includes sex-specific phenos)
all.out.0.nocov <- scanone(sfon, method="hk", pheno.col = c(names(ph.no.cov), ph.sex.sp))

# Genome-wide permutations
all.out.0.nocov.perm <- scanone(sfon, method="hk", pheno.col= c(names(ph.no.cov), ph.sex.sp),
                                n.perm=num.perm, verbose=T,
                                n.cluster = num.cluster)


####### 2B SINGLE QTL, consider covariate #####
# scanone (consider covariate)
sex <- as.numeric(pull.pheno(sfon, "sex") == "M") #create numeric sex variable (fem 0 ; male 1)
out.am <- scanone(sfon, method="hk", addcov=sex, 
                  pheno.col = c(names(ph.yes.cov),names(ph.no.cov))
) # additive model
out.im <- scanone(sfon, method="hk", addcov=sex, 
                  intcovar=sex, 
                  pheno.col = c(names(ph.yes.cov),names(ph.no.cov))
) # 'full' model

# Genome-Wide permutations
# set seed to directly compare full model and full-additive model
set.seed(54955149)
operm.am <- scanone(sfon, addcovar=sex, n.perm=num.perm,
                    pheno.col = c(names(ph.yes.cov),names(ph.no.cov)),
                    n.cluster = num.cluster)
set.seed(54955149)
operm.im <- scanone(sfon, addcovar=sex, intcovar=sex, n.perm=num.perm,
                    pheno.col = c(names(ph.yes.cov),names(ph.no.cov)),
                    n.cluster = num.cluster)
#note: at ~ 100perms/2h; therefore to scale make sure to use n.cluster


## Chromosome-wide significance using for loop for pheno and chr, 
## includes if/else statement to use the covariate of sex when required
## collects scanone objects and the chromosome-level p = (chr.wide.pval) sig LOD

# select your phenotypes and chromosomes
all.phenos <- c(names(ph.no.cov), ph.sex.sp, names(ph.yes.cov))
selected.chrs <- 1:42

# set your pval
chr.wide.pval <- 0.01

# set NULLs
chr.scan = NULL # a scanone object rewritten each chromosome
chr.scan.perm = NULL # scanone object with permutations rewritten each chromosome
chr.sig = NULL # this is the LOD level for each chromosome at p = 0.05
scanone.mods <- list() # this is a list of scanone objects, one for each pheno x chr
pheno.names = NULL # this is a vector of the names of tested phenotypes
pheno.sig.lod.per.chr = NULL # this is the LOD sig level for each pheno x chr

# Loop across phenos, and chromosomes
for(pheno in all.phenos) {
  chr.sig = NULL # set NULL each pheno loop
  for(chr in selected.chrs) {
if (pheno %in% names(ph.yes.cov)) {
  # use covariate if needed
  chr.scan <- scanone(sfon, method = "hk", pheno.col=pheno, chr=chr, addcov=sex)
  chr.scan.perm = scanone(sfon, method="hk", 
                          pheno.col=pheno, chr=chr, addcov=sex,
                          n.perm=num.perm, n.cluster = num.cluster 
                          , verbose=T)
  chr.sig = c(chr.sig, summary(chr.scan.perm, chr.wide.pval))
  scanone.mods[[paste(pheno, "_", chr, sep="")]] <- chr.scan
} else {
  # don't use covariate
  chr.scan <- scanone(sfon, method = "hk", pheno.col=pheno, chr=chr)
  chr.scan.perm = scanone(sfon, method="hk", 
                          pheno.col=pheno, chr=chr, 
                          n.perm=num.perm, n.cluster = num.cluster 
                          , verbose=T)
  chr.sig = c(chr.sig, summary(chr.scan.perm, chr.wide.pval))
  scanone.mods[[paste(pheno, "_", chr, sep="")]] <- chr.scan
}
}
pheno.names = c(pheno.names, names(sfon$pheno[pheno])) # obtain the name of the tested pheno in this loop
pheno.sig.lod.per.chr = cbind(pheno.sig.lod.per.chr, chr.sig)
}
colnames(pheno.sig.lod.per.chr) <- pheno.names #gives names to the average significance levels
pheno.sig.lod.per.chr
scanone.mods


str(scanone.mods)


#####2C SEX AS A BINARY TRAIT####
sfon$pheno <- cbind(sfon$pheno, binary=sex)
out.bin <- scanone(sfon, pheno.col="binary",
                   model="binary")

# Genome-wide significance
operm.bin <- scanone(sfon, pheno.col="binary", 
                     model="binary", n.perm=num.perm, n.cluster = num.cluster)


#####3D INTERVAL ESTIMATES#####
#?bayesint # uses output from scanone
#use in C Sauvage's analysis: 
# bayesint(outall.hk, 5, 0.95, expandtomarkers=T) # chr = 5, prob=0.95, allow interval to flanking markers
# out.boot = scanoneboot(raw.data, chr = 1:40, n.boot = 100) #cross = 'raw.data'
#?scanoneboot # get estimated confidence interval for the location of a QTL

#Try 1.5 LOD interval

# Histogram of the estimated QTL locations in 1000 bootstrap replicates
#hist(out.boot, main="Bootstrap values of QTL position", xlab="Estimated QTL location (cM)",
#     ylab="", yaxt="n")

#rug(raw.data$geno[[5]]$map) # not sure what this is showing,


### SAVE OUT PERMUTATION OBJECT ####
#save.image(file = "sfon_02_output.RData") # save out existing data so that it can be reloaded without running all again

#END