# Analysis of the phenotypic QTL from Sauvage et al with the new map from Sutherland et al.
# analysis v2.0
# Part 02/03 (permutations on cluster computer)
# in R/QTL (cp cross type)

# rm(list=ls())

# install libraries
require(qtl)
# install.packages("rlecuyer")
require(rlecuyer) #random number generator in parallel
require(plyr)

# set working directory for the cluster computer
setwd("~/Documents/sfqtl_macpro_2016-09-16/sfon_pqtl")

# Load Part 1 results:
load("02_data/sfon_01_output.RData")

## Also can load subset Part 1 results
#load("02_data/sfon_01_output_subset_only_nnxnp_and_mf_sep.RData") # subset markers
load("02_data/sfon_01_output_subset_only_nnxnp_and_hkxhk.RData") # subset markers
# sfon <- sfon.males
# sfon <- sfon.females
# sfon <- sfon_only_nnxnp

# sfon <- sfon_limited

## note will need to rerun
# sfon <- est.rf(sfon)
# sfon <- calc.genoprob(sfon, step = 0, error.prob =0.0001) #hmm to calc prob of true underlying genos given data
sfon

# Set permutation variables (perm should be 1000):
num.perm = 10000
num.cluster = 20

sex <- as.numeric(pull.pheno(sfon, "sex") == "M") #create numeric sex variable (fem 0 ; male 1)

# Light analysis (with sex info)
out <- scanone(sfon, method="hk", pheno.col = c(names(ph.yes.cov),names(ph.no.cov))
               , addcov=sex)
set.seed(54955149)
operm.out <- scanone(sfon, n.perm=num.perm,
                     pheno.col = c(names(ph.yes.cov),names(ph.no.cov)),
                     n.cluster = num.cluster
                     , addcovar=sex)


### SEX SPECIFIC ANALYSIS (TEMP) ###
# scanone (consider covariate)
out <- scanone(sfon, method="hk", pheno.col = c(names(ph.yes.cov),names(ph.no.cov)))
set.seed(54955149)
operm.out <- scanone(sfon, n.perm=num.perm,
                    pheno.col = c(names(ph.yes.cov),names(ph.no.cov)),
                    n.cluster = num.cluster)

names(out)[3:29]
summary(out, perms=operm.out, format="tabByCol", alpha=0.15, pvalues=T)

# start by setting your variables
POI <- "" #set your pheno
scanone.mod <- out
scanone.perms <- operm.out
lod.col <- which(names(scanone.mod) == POI) - 2 #-2 is to account for chr and pos

# then plot
plot(scanone.mod, alternate.chrid=T, lodcolumn=lod.col,
     ylim = c(0,10), 
     ylab = POI,
     xlab = "",
     bandcol="gray70")
abline(h=summary(scanone.perms[,POI], 0.05), lty=1)

### END SEX SPECIFIC ANALYSIS (TEMP) ###



# Create sex variable
sex <- as.numeric(pull.pheno(sfon, "sex") == "M") #create numeric sex variable (fem 0 ; male 1)

####### 2A Single QTL, no covariate #####
# Create scanone object for phenos not requiring covariate (includes sex-specific phenos)
all.out.0.nocov <- scanone(sfon, method="hk", pheno.col = c(names(ph.no.cov), ph.sex.sp))

# Genome-wide permutations
all.out.0.nocov.perm <- scanone(sfon, method="hk", pheno.col= c(names(ph.no.cov), ph.sex.sp),
                                n.perm=num.perm, verbose=T,
                                n.cluster = num.cluster)


####### 2B SINGLE QTL, consider covariate #####
# scanone (consider covariate)
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


## 2C Chromosome-wide scanone and permutations ####
## includes if/else statement to use the covariate of sex when required
# select your phenotypes and chromosomes
all.phenos <- c(names(ph.no.cov), ph.sex.sp, names(ph.yes.cov))
selected.chrs <- 1:42

# set your pval
chr.wide.pval <- 0.01

# find which phenotypes are numeric
nums <- sapply(sfon$pheno, is.numeric) # find out which phenotypes to assess (are numeric)

print(paste(c("Using", num.perm, "permutations")), quote=F)

results <- capture.output(
  for(m in which(nums)) {
    x.hk <- NULL
    x.perm <- NULL
    i=NULL
    print(paste(c(m, names(sfon$pheno[m])))) # print the phenotype being tested
    if(names(nums[m]) %in% c(names(ph.no.cov), ph.sex.sp)) {
      print("Not incl covariate") # do not use covariate
      for(i in 1:nchr(sfon)){
        x.hk=scanone(sfon, method="hk", pheno.col=m, chr=i)
        x.perm=scanone(sfon, method="hk", n.perm=num.perm, n.cluster=num.cluster, pheno.col=m, chr=i, verbose=F)
        print(paste(c("Testing chromosome", i))) # print the chromosome being tested 
        print(summary(x.hk, perms=x.perm, alpha=0.01, pvalues=T, format="tabByCol")) #find the LOD peak for this
      }
    }
    else if(names(nums[m]) %in% names(ph.yes.cov)){
      print("Including covariate")
      for(i in 1:nchr(sfon)){
        x.hk=scanone(sfon, method="hk", pheno.col=m, chr=i, addcov=sex)
        x.perm=scanone(sfon, method="hk", n.perm=num.perm, n.cluster=num.cluster, pheno.col=m, chr=i, addcov=sex, verbose=F)
        print(paste(c("Testing chromosome", i))) # print the chromosome being tested 
        print(summary(x.hk, perms=x.perm, alpha=0.01, pvalues=T, format="tabByCol")) #find the LOD peak for this
      }
    }
  }
)


write.csv(results, file = "new_chr_wide_p0.01_w_10000perms.csv", quote=F, row.names=F)

#####2D SEX AS A BINARY TRAIT####
sfon$pheno <- cbind(sfon$pheno, binary=sex)
out.bin <- scanone(sfon, pheno.col="binary",
                   model="binary")

# Genome-wide significance
operm.bin <- scanone(sfon, pheno.col="binary", 
                     model="binary", n.perm=num.perm, n.cluster = num.cluster)

### SAVE OUT PERMUTATION OBJECT ####
#save.image(file = "sfon_02_output.RData") # save out existing data so that it can be reloaded without running all again
## save.image(file = "02_data/sfon_02_output_gw_1000perms_nnxnp_only.RData")
#END