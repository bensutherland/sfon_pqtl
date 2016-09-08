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

##### is this package necessary ?####
#install.packages("snow")
#library(snow)

# set working directory for the cluster computer

# Load Part 1 results:
load("sfon_01_output.RData")

# Set permutation variables (perm should be 1000):
num.perm = 1000
num.cluster = 20

####### 2A Single QTL, no covariate #####
# scanone (no covar), includes sex-specific phenos
all.out.0.nocov <- scanone(sfon, method="hk", pheno.col = c(names(ph.no.cov), ph.sex.sp))
# Genome-Wide permutations
all.out.0.nocov.perm <- scanone(sfon, method="hk", pheno.col= c(names(ph.no.cov), ph.sex.sp),
                                n.perm=num.perm, verbose=T,
                                n.cluster = num.cluster)



#######2B SINGLE QTL, consider covariate#####
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
save.image(file = "sfon_02_output.RData") # save out existing data so that it can be reloaded without running all again

#END