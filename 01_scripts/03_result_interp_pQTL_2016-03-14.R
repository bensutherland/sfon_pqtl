# Analysis of the phenotypic QTL from Sauvage et al with the new map from Sutherland et al.
# analysis v2.0
# Part 03/03 (result interpretation and figure generation)
# in R/QTL (cp cross type)

# rm(list=ls())

# install libraries
require(qtl)

# Load data from permutation tests
load("02_data/sfon_02_output_gw_and_chr.RData") # both chromosome and genome-wide significance THIS
#or
load("02_data/sfon_02_output.RData") # currently set up below USE THIS ONE

# other datasets: (to be cleaned up)
#load("02_data/sfon_01_output.RData")
#load("02_data/sfon_02_output_chromosome-wide_p0.01_1000perms.RData") # 
#load("02_data/sfon_02_output_chromosome-wide_p0.01_10000perms.RData") # 
#load("02_data/sfon_02_output_gw_1000perms_nnxnp_only.RData") # only female-specific markers

# For most of the plotting (following) we need to simulate genotypes given observed marker data
sfon <- sim.geno(sfon, step=2.5,
                  error.prob=0.001, n.draws=256) #note that 4-way cross w sex-specific map assumes constant ratio of female:male recomb rates within inter-marker intervals
str((sfon$geno[[1]])) # note that there is now probabilities and draws added to each LG

# Data overview
# these are your phenotypes
names(ph.no.cov)
ph.sex.sp
names(ph.yes.cov)

# Need to setup sex as a dataframe for PVE estimates
sex.df <- as.data.frame(sex)


######### 3A RESULT INTERP, SINGLE QTL, no COV #######
# Note: uses perms from part 2
names(ph.no.cov)
summary(all.out.0.nocov, perms=all.out.0.nocov.perm, alpha=0.15, format="tabByCol", pvalues=TRUE)


## Plot your LOD curves and genome-wide significance
# Start by setting your variables
POI <- "ghr" #set your pheno
scanone.mod <- all.out.0.nocov
scanone.perms <- all.out.0.nocov.perm

# Then run the following
lod.col <- which(names(scanone.mod) == POI) - 2 #-2 is to account for chr and pos
plot(scanone.mod, alternate.chrid=T, lodcolumn=lod.col,
     ylim = c(0,8), 
     ylab = POI,
     xlab = "",
     bandcol="gray70")
abline(h=summary(scanone.perms[,POI], 0.05), lty=1)

# If saving, scanone objects as 11 x 5


######### 3B RESULT INTERP, SINGLE QTL, consider COV #######
names(ph.yes.cov)

# BACKGROUND INFO  & DATA SETUP #
head(operm.am) #Note: no pvals assoc. w scanone obj, pvals intro'd by 'summary'
out.im[1,] # each model obj. contains LOD scores, which can be compared
out.am[1,]
out.im[1,] - out.am[1,]
## for easy use, combine full and interactive models
out.fmim <- c(out.im, out.im - out.am, 
              labels=c("f","i"))  # error expected?
names(out.fmim) # note the positions of the full and interactive values

# for ease, combine full and interactive models with the permutations applied
operm.fmim <- cbind(operm.im, operm.im - operm.am,
                    labels=c("f","i"))
colnames(operm.im) # NOTE that this contains all phenotypes, including those that do not require sex as cov.
colnames(operm.fmim)

# COMPARE RESULTS TO IDENTIFY SIG INTERACTION #
summary(operm.fmim, alpha=0.05) # provides the LOD threshold for p ≤ 0.05

# Check for significance in the full model and in the interactive model alone
summary(out.fmim, perms=operm.fmim, format="tabByCol", alpha=0.15, pvalues=TRUE)
# If high LOD comes from the interaction, rather than additive portion, sig ifx.
# If trait has no sig ifx, check sig in additive model only (interaction uses up power)
summary(out.am, perms=operm.am, format="tabByCol", alpha=0.1, pvalues=T)

# PLOT TO INSPECT INTERACTION #
# View the LOD from the full, interactive and additive models
# for example, a significant interaction was found with osmo.delta on chr13
names(out.fmim)
trait.lod = 9 #choose the trait lod column to plot
plot(out.im, out.im-out.am, out.am, ylab = "LOD score",
     col=c("blue", "red", "black"), alternate.chrid=T, lodcolumn=trait.lod,
     #main = "osmo.delta"
     )
legend(x = "topright", legend=c("full","interact","additive"), 
       #col = c("blue", "red", "black"), 
       fill=c("blue", "red", "black"))


## PLOT SCANONE ADDITIVE MODELS ##
#find your LOD column (start counting after chr and pos)
names(out.am) # needs covariate

## REGULAR
# start by setting your variables
POI <- "weight_liver.g" #set your pheno
scanone.mod <- out.am
scanone.perms <- operm.am
lod.col <- which(names(scanone.mod) == POI) - 2 #-2 is to account for chr and pos

# then plot
plot(scanone.mod, alternate.chrid=T, lodcolumn=lod.col,
     ylim = c(0,10), 
     ylab = POI,
     xlab = "",
     bandcol="gray70")
abline(h=summary(scanone.perms[,POI], 0.05), lty=1)


#Needs some special playing around to get the full model version rather than the additive model version
names(out.fmim)
plot(out.fmim, alternate.chrid = T, lodcolumn="weight_liver.g.f")



## SPECIAL TWO TIME POINT PLOT ##
# start by setting your variables
POI <- "weight.g_0509" #set your pheno
POINAME <- "weight (g)"
scanone.mod <- out.am
scanone.perms <- operm.am
lod.col <- which(names(scanone.mod) == POI) - 2 #-2 is to account for chr and pos

# then plot
plot(scanone.mod, alternate.chrid=T, lodcolumn=lod.col,
     ylim = c(0,6), 
     ylab = POINAME,
     xlab = "",
     bandcol="gray70"
     #, add = T
     , col = "red"
     , lty = 2
     )
abline(h=summary(scanone.perms[,POI], 0.05), lty=3)

# start by setting your variables
POI <- "weight.g_1109" #set your pheno
scanone.mod <- out.am
scanone.perms <- operm.am
lod.col <- which(names(scanone.mod) == POI) - 2 #-2 is to account for chr and pos

# then plot
plot(scanone.mod, alternate.chrid=T, lodcolumn=lod.col,
     ylim = c(0,6), 
     ylab = POINAME,
     xlab = "",
     bandcol="gray70"
     , add = T
     #, col = "blue"
     , lty = 1
)
abline(h=summary(scanone.perms[,POI], 0.05), lty=3)

legend("topright", legend = c("T1","T3"), fill = c("red", "black"))

### END SPECIAL TWO PLOT ###


######## 3C EXTRACT RESULTS OF CHR-WIDE SIG ######
summary(scanone.mods[["weight.g_0509_1"]], threshold = 1.2) # for example
summary(scanone.mods[["weight.g_0509_1"]], threshold = pheno.sig.lod.per.chr[1,"weight.g_0509"])

names(scanone.mods)
summary(scanone.mods[["cort.delta_6"]], threshold = pheno.sig.lod.per.chr[1,"cort.delta"])


# Identify and extract chromosome-wide significant QTL from scanone object:

chr.wide.output <- NULL

chr.wide.output <- 
  capture.output(
  for(pheno in all.phenos) {
    for(chr in selected.chrs) {
      paste(
        cat(c(pheno, chr, "result"), sep = "_"),
        print(summary(scanone.mods[[paste(pheno, "_", chr, sep="")]], threshold = pheno.sig.lod.per.chr[chr, pheno]))
        )
      }
  }
)

str(chr.wide.output)
tail(chr.wide.output, n = 5)

write.csv(chr.wide.output, file = "02_data/chr_wide_p0.01_w_10000perms.csv", quote = F, row.names=F)
# then in terminal:
# grep -A1 'lod' 02_data/chr_wide_p0.01_w_10000perms.csv | grep -v '--' - | less 
# or redirect to file

######### 3D RESULT INTERP, SEX AS A BINARY TRAIT #######
# Bare bones plot
plot(out.bin, col="red", ylab="LOD score",
     alternate.chrid=TRUE)

# Better plot
plot(out.bin, col="red", ylab="LOD score",
     alternate.chrid=TRUE
     , bandcol = "lightgrey"
    )
# add significance line
abline(h = summary(operm.bin, alpha=0.05), lty = 2)

# here are the significance values
summary(operm.bin, alpha=0.05) # provides the LOD threshold for p ≤ 0.05
summary(out.bin, perms=operm.bin, alpha=0.05, pvalues=TRUE) # provides the highest position per chromosome where LOD > thresh

# clearly sex-linked chromosome at 35
# save out as 10 x 4

######### 3E PLOT DATA, sex markers #######
# find marker with highest LOD
find.marker(cross = sfon, chr = 35, pos = 109) 
# Note: if the marker is not fully informative, genotype plotting will vary each time due to imputation of genotypes
# so it may be more informative to find the closest fully informative marker as long as the lod is > thresh

plot.pxg(sfon, marker="43393", pheno.col = "binary") #infer=F doesn't prevent imputation         
# infer = T means that missing genos are filled with a random imputation and plotted in red

effectplot(sfon, pheno.col= "binary", mname1="83995") # still uses imputation


######### 3F PLOT DATA, EFFECT PLOTS #######

# This section needs some work #

# find a trait and marker with high lod, with an additive effect
summary(out.am, perms=operm.am, format="tabByCol", alpha=0.1, pvalues=T)
summary(out.fmim, perms=operm.fmim, format="tabByCol", alpha=0.1, pvalues=T) # make sure no interact
# e.g. for weight.g_0509 , consider chr 5 marker 85980 

par(mfrow=c(1,2), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

#example with covariate, but using the actual data, not imputed data
# first lets find the actual data
sfon$geno[[4]]$data[1:5,1:5] # gives the genotype codes from the first five markers on chr 5 for the first five indiv (out of 170)
#Important## if using mark1 make sure to use geno1 or else arbitrary rqtl genotype labels will be used (i.e. first 2)

# as an example, check out a fully informative marker (with four distinct genotypes):
par(mfrow=c(1,1))
marker.count <- which(names(sfon$geno[[4]]$data[1,]) == "7187") #find a marker by name on chr 4
effectplot(sfon, mname1="7187",  
           #mname2="Sex",
           mark1=sfon$geno[[4]]$data[,marker.count], 
           #geno1=c("AG,AG","GA,AG","AG,AA","GA,AA"), #appropriate genos are added from stacks.markers.info.xlsx
           #mark2=sex, geno2=c("F","M"), 
           pheno.col = "weight.g_0709", 
           #add.legend=F, 
           ylim = c(100,300),
           ylab = "weight at T2",
           main = "",
           xlab = ""
) 
#use geno1 to give the actual genotypes

# is there a significant chr-wide interaction


# automatic capture of the genotypes using 'mname' without assigning mark1 will
# use the genotype data that has been imputed, and thus may be slightly dependent on
# the sim.geno run.
effectplot(sfon, mname1 = "7187", mname2 = "Sex",
           mark2=sex, geno2 = c("F","M"),
           pheno.col = "weight.g_0709",
           add.legend=T,
           ylim = c(100,300),
           main = "")

########## 3G PVE Genome Wide ######
# for this use the following formulae (per trait)
?makeqtl()
qtl.trait <- makeqtl(sfon, chr=c(), pos=c())
?fitqtl()
out.fq <- fitqtl(sfon, qtl=qtl.trait, formula=y~Q1+Q2+Q3) # put all of the QTL in one formula
summary(out.fq) # this will give the PVE
plot(qtl.object) # this will put the QTL on your genetic map

# example with weight.g_0509

qtl.weight.g_0509 <- makeqtl(sfon, chr = c("5","20"), pos = c(261,162)
                             , what="draws" # to use prob, requires imputation using hk
)
summary(qtl.weight.g_0509)
out.weight.g_0509.fq <- fitqtl(sfon, pheno.col="weight.g_0509", qtl=qtl.weight.g_0509, formula=y~Q1+Q2
                               #, covar=sex.df
)

summary(out.weight.g_0509.fq) # need to put in the covariate? doesn't seem to affect anything

# example with weight.g_0709
qtl.weight.g_0709 <- makeqtl(sfon, chr=c("4","5"), pos = c(28.3, 198.5), what = "draws")
out.weight.g_0709 <- fitqtl(sfon, pheno.col="weight.g_0709", qtl=qtl.weight.g_0709, formula=y~Q1+Q2)
summary(out.weight.g_0709)

# as above, but including the sex covariate (additive)
sex.df <- as.data.frame(sex)
qtl.weight.g_0709_sex <- makeqtl(sfon, chr=c("4","5"), pos = c(28.3, 198.5), what = "draws")
out.weight.g_0709_sex <- fitqtl(sfon, pheno.col="weight.g_0709", qtl=qtl.weight.g_0709_sex
                                , formula=y~Q1+Q2+sex, covar=sex.df)
summary(out.weight.g_0709_sex)




# example to find pve from a single qtl
qtl.weight.g_1109 <- makeqtl(sfon, chr=c("20"), pos = c(169), what = "draws")
out.weight.g_1109 <- fitqtl(sfon, pheno.col="weight.g_1109", qtl=qtl.weight.g_1109, formula=y~Q1)
summary(out.weight.g_1109)

# TCS_T1.T2
qtl.TCS_T1.T2 <- makeqtl(sfon, chr=c("9"), pos = c(184), what = "draws")
out.TCS_T1.T2 <- fitqtl(sfon, pheno.col="TCS_T1.T2", qtl=qtl.TCS_T1.T2, formula=y~Q1)
summary(out.TCS_T1.T2)

# TCS_T2.T3
qtl.TCS_T2.T3 <- makeqtl(sfon, chr=c("3","36"), pos = c(35.7, 31.4), what = "draws")
out.TCS_T2.T3 <- fitqtl(sfon, pheno.col="TCS_T2.T3", qtl=qtl.TCS_T2.T3, formula=y~Q1+Q2)
summary(out.TCS_T2.T3)

# TCS_T1.T3
qtl.TCS_T1.T3 <- makeqtl(sfon, chr=c("36"), pos = c(61.8), what = "draws")
out.TCS_T1.T3 <- fitqtl(sfon, pheno.col="TCS_T1.T3", qtl=qtl.TCS_T1.T3, formula=y~Q1)
summary(out.TCS_T1.T3)

# leng.cm_0509
qtl.leng.cm_0509 <- makeqtl(sfon, chr=c("5","34"), pos = c(261,110), what = "draws")
out.leng.cm_0509 <- fitqtl(sfon, pheno.col="leng.cm_0509", qtl=qtl.leng.cm_0509, formula=y~Q1+Q2)
summary(out.leng.cm_0509)

# leng.cm_0709
qtl.leng.cm_0709 <- makeqtl(sfon, chr=c("4","5"), pos = c(115,261), what = "draws")
out.leng.cm_0709 <- fitqtl(sfon, pheno.col="leng.cm_0709", qtl=qtl.leng.cm_0709, formula=y~Q1+Q2)
summary(out.leng.cm_0709)

# leng.cm_1109
qtl.leng.cm_1109 <- makeqtl(sfon, chr=c("20"), pos = c(169), what = "draws")
out.leng.cm_1109 <- fitqtl(sfon, pheno.col="leng.cm_1109", qtl=qtl.leng.cm_1109, formula=y~Q1)
summary(out.leng.cm_1109)

# condit.fact_T3
qtl.condit.fact_T3 <- makeqtl(sfon, chr=c("20"), pos = c(150), what = "draws")
out.condit.fact_T3 <- fitqtl(sfon, pheno.col="condit.fact_T3", qtl=qtl.condit.fact_T3, formula=y~Q1)
summary(out.condit.fact_T3)

# weight_liver.g
qtl.weight_liver.g <- makeqtl(sfon, chr=c("35"), pos = c(24.5), what = "draws")
out.weight_liver.g <- fitqtl(sfon, pheno.col="weight_liver.g", qtl=qtl.weight_liver.g, formula=y~Q1)
summary(out.weight_liver.g)

# osmo.delta
qtl.osmo.delta <- makeqtl(sfon, chr=c("13"), pos = c(26.5), what = "draws")
out.osmo.delta <- fitqtl(sfon, pheno.col="osmo.delta", qtl=qtl.osmo.delta, formula=y~Q1)
summary(out.osmo.delta)

# fem.egg.diam
qtl.fem.egg.diam <- makeqtl(sfon, chr=c("5"), pos = c(66.8), what = "draws")
out.fem.egg.diam <- fitqtl(sfon, pheno.col="fem.egg.diam", qtl=qtl.fem.egg.diam, formula=y~Q1)
summary(out.fem.egg.diam)

# plasma.chlor
qtl.plasma.chlor <- makeqtl(sfon, chr=c("31"), pos = c(113.9), what = "draws")
out.plasma.chlor <- fitqtl(sfon, pheno.col="plasma.chlor", qtl=qtl.plasma.chlor, formula=y~Q1)
summary(out.plasma.chlor)

# ghr
qtl.ghr <- makeqtl(sfon, chr=c("24"), pos = c(138), what = "draws")
out.ghr <- fitqtl(sfon, pheno.col="ghr", qtl=qtl.ghr, formula=y~Q1)
summary(out.ghr)

##### 3H PVE chromosome-wide #######

## Phenotypes not needing covariate sex ##
# TCS_T1.T2
qtl.chr.TCS_T1.T2 <- makeqtl(sfon, chr=c("9"), pos = c(184), what = "draws")
out.chr.TCS_T1.T2 <- fitqtl(sfon, pheno.col="TCS_T1.T2", qtl=qtl.chr.TCS_T1.T2, formula=y~Q1
                            , get.ests = TRUE
                            )
summary(out.chr.TCS_T1.T2)

# TCS_T2.T3
qtl.chr.TCS_T2.T3 <- makeqtl(sfon, chr=c("3","9","36"), pos = c(35.7, 139, 31.4), what = "draws")
out.chr.TCS_T2.T3 <- fitqtl(sfon, pheno.col="TCS_T2.T3", qtl=qtl.chr.TCS_T2.T3, formula=y~Q1+Q2+Q3
                            , get.ests = TRUE
                            )
summary(out.chr.TCS_T2.T3)

# TCS_T1.T3
qtl.chr.TCS_T1.T3 <- makeqtl(sfon, chr=c("36"), pos = c(61.8), what = "draws")
out.chr.TCS_T1.T3 <- fitqtl(sfon, pheno.col="TCS_T1.T3", qtl=qtl.chr.TCS_T1.T3, formula=y~Q1)
summary(out.chr.TCS_T1.T3)

# condit.fact_T1
qtl.chr.condit.fact_T1 <- makeqtl(sfon, chr=c("16"), pos = c(89.9), what = "draws")
out.chr.condit.fact_T1 <- fitqtl(sfon, pheno.col="condit.fact_T1", qtl=qtl.chr.condit.fact_T1, formula=y~Q1)
summary(out.chr.condit.fact_T1)

# condit.fact_T2
qtl.chr.condit.fact_T2 <- makeqtl(sfon, chr=c("39"), pos = c(46.9), what = "draws")
out.chr.condit.fact_T2 <- fitqtl(sfon, pheno.col="condit.fact_T2", qtl=qtl.chr.condit.fact_T2, formula=y~Q1)
summary(out.chr.condit.fact_T2)

# condit.fact_T3
qtl.chr.condit.fact_T3 <- makeqtl(sfon, chr=c("20"), pos = c(150), what = "draws")
out.chr.condit.fact_T3 <- fitqtl(sfon, pheno.col="condit.fact_T3", qtl=qtl.chr.condit.fact_T3, formula=y~Q1)
summary(out.chr.condit.fact_T3)

# chlor.delta
qtl.chr.chlor.delta <- makeqtl(sfon, chr=c("38"), pos = c(43.3), what = "draws")
out.chr.chlor.delta <- fitqtl(sfon, pheno.col="chlor.delta", qtl=qtl.chr.chlor.delta, formula=y~Q1)
summary(out.chr.chlor.delta)

# fem.egg.diam
qtl.chr.fem.egg.diam <- makeqtl(sfon, chr=c("5"), pos = c(66.8), what = "draws")
out.chr.fem.egg.diam <- fitqtl(sfon, pheno.col="fem.egg.diam", qtl=qtl.chr.fem.egg.diam, formula=y~Q1)
summary(out.chr.fem.egg.diam)

# male.sperm.diam
qtl.chr.male.sperm.diam <- makeqtl(sfon, chr=c("24"), pos = c(47.8), what = "draws")
out.chr.male.sperm.diam <- fitqtl(sfon, pheno.col="male.sperm.diam", qtl=qtl.chr.male.sperm.diam, formula=y~Q1)
summary(out.chr.male.sperm.diam)

# ghr
qtl.chr.ghr <- makeqtl(sfon, chr=c("24"), pos = c(138), what = "draws")
out.chr.ghr <- fitqtl(sfon, pheno.col="ghr", qtl=qtl.chr.ghr, formula=y~Q1)
summary(out.chr.ghr)

# # Phenos needing sex as covariate
# example how to include sex as additive covariate
# qtl.pheno <- makeqtl(sfon, chr=c(), pos = c(), what = "draws")
# out.pheno <- fitqtl(sfon, pheno.col="", qtl=qtl.pheno
#                                , formula=y~Q1+Q2+sex, covar=sex.df)
# summary(out.pheno)

# weight.g_0509
qtl.chr.weight.g_0509 <- makeqtl(sfon, chr=c("4","5","20"), pos = c(28.3,261,162), what = "draws")
out.chr.weight.g_0509 <- fitqtl(sfon, pheno.col="weight.g_0509", qtl=qtl.chr.weight.g_0509, formula=y~Q1+Q2+Q3+sex, covar=sex.df)
summary(out.chr.weight.g_0509)

# weight.g_0709
qtl.chr.weight.g_0709 <- makeqtl(sfon, chr=c("4","5","20"), pos = c(28.3,198,162), what = "draws")
out.chr.weight.g_0709 <- fitqtl(sfon, pheno.col="weight.g_0709", qtl=qtl.chr.weight.g_0709, formula=y~Q1+Q2+Q3+sex, covar=sex.df)
summary(out.chr.weight.g_0709)

# weight.g_1109
qtl.chr.weight.g_1109 <- makeqtl(sfon, chr=c("20","24"), pos = c(169,79.8), what = "draws")
out.chr.weight.g_1109 <- fitqtl(sfon, pheno.col="weight.g_1109", qtl=qtl.chr.weight.g_1109, formula=y~Q1+Q2+sex, covar=sex.df)
summary(out.chr.weight.g_1109)

# leng.cm_0509
qtl.chr.leng.cm_0509 <- makeqtl(sfon, chr=c("5","20","27","34"), pos = c(261,162,45.1,110), what = "draws")
out.chr.leng.cm_0509 <- fitqtl(sfon, pheno.col="leng.cm_0509", qtl=qtl.chr.leng.cm_0509, formula=y~Q1+Q2+Q3+Q4+sex, covar=sex.df)
summary(out.chr.leng.cm_0509)

# leng.cm_0709
qtl.chr.leng.cm_0709 <- makeqtl(sfon, chr=c("3","4","5","20","34"), pos = c(83,115,261,169,112), what = "draws")
out.chr.leng.cm_0709 <- fitqtl(sfon, pheno.col="leng.cm_0709", qtl=qtl.chr.leng.cm_0709, formula=y~Q1+Q2+Q3+Q4+Q5+sex, covar=sex.df)
summary(out.chr.leng.cm_0709)

# leng.cm_1109
qtl.chr.leng.cm_1109 <- makeqtl(sfon, chr=c("3","20","24"), pos = c(52.9,169,79.8), what = "draws")
out.chr.leng.cm_1109 <- fitqtl(sfon, pheno.col="leng.cm_1109", qtl=qtl.chr.leng.cm_1109, formula=y~Q1+Q2+Q3+sex, covar=sex.df)
summary(out.chr.leng.cm_1109)

# weight_liver.g
qtl.chr.weight_liver.g <- makeqtl(sfon, chr=c("1","2","35"), 
                                  pos = c(148,318,46.5), what = "draws")
out.chr.weight_liver.g <- fitqtl(sfon, pheno.col="weight_liver.g", qtl=qtl.chr.weight_liver.g, formula=y~Q1+Q2+Q3)
summary(out.chr.weight_liver.g)

# hematocr
qtl.chr.hematocr <- makeqtl(sfon, chr=c("4","25"), pos = c(22.6,139), what = "draws")
out.chr.hematocr <- fitqtl(sfon, pheno.col="hematocr", qtl=qtl.chr.hematocr, formula=y~Q1+Q2+sex, covar=sex.df)
summary(out.chr.hematocr)

# cort.delta
qtl.chr.cort.delta <- makeqtl(sfon, chr=c("6"), pos = c(114), what = "draws")
out.chr.cort.delta <- fitqtl(sfon, pheno.col="cort.delta", qtl=qtl.chr.cort.delta, formula=y~Q1+sex, covar=sex.df)
summary(out.chr.cort.delta)

# osmo.delta
qtl.chr.osmo.delta <- makeqtl(sfon, chr=c("40"), pos = c(60.3), what = "draws")
out.chr.osmo.delta <- fitqtl(sfon, pheno.col="osmo.delta", qtl=qtl.chr.osmo.delta, formula=y~Q1+sex, covar=sex.df)
summary(out.chr.osmo.delta)


#### 3I Effect sizes per mtype per sig marker/trait #####
# import dataframe of important markers and chromosome
# Collect the columns: phenotype; chr; and mname -- then put into a .csv
sig_mname_chr_pheno.df <- as.data.frame(read.csv(file = "02_data/sig_pheno_chr_mname_from_excel.csv", header = T))
head(sig_mname_chr_pheno.df)
sig_mname_chr_pheno.df[] <- lapply(sig_mname_chr_pheno.df, as.character) #convert the dataframe cols to char
str(sig_mname_chr_pheno.df)

# Use aggregate to calculate the mean per marker type per sig marker/trait combos
phenotype=NULL
mname=NULL
chr=NULL
for(i in 1:length(sig_mname_chr_pheno.df[,1])) {
  mname=sig_mname_chr_pheno.df[i,"mname"]
  chr=sig_mname_chr_pheno.df[i,"chr"]
  phenotype=sig_mname_chr_pheno.df[i,"phenotype"]
  print(c(mname,chr,phenotype))
  print(aggregate(sfon$pheno[phenotype], by=list(sfon$geno[[chr]]$data[,mname]), FUN=mean, na.rm=T))
}

# in case want to test out, note that the na.rm=T is critical here
# it is not clear why na.rm is critical, because sometimes works when there are NA values and sometimes doesn't
aggregate(sfon$pheno$TCS_T1.T2, by=list(sfon$geno[[16]]$data[,"118085"]), FUN=mean
          , na.rm=T
          )

# as above, but also including sex in model when needed
phenotype=NULL
mname=NULL
chr=NULL
for(i in 1:length(sig_mname_chr_pheno.df[,1])) {
  mname=sig_mname_chr_pheno.df[i,"mname"]
  chr=sig_mname_chr_pheno.df[i,"chr"]
  phenotype=sig_mname_chr_pheno.df[i,"phenotype"]
  print(c(mname,chr,phenotype))
  if (phenotype %in% names(ph.yes.cov)) {
    print("sex included")
    print(aggregate(sfon$pheno[phenotype], by=list(sfon$geno[[chr]]$data[,mname], sfon$pheno$sex), FUN=mean, na.rm=T))
  } else if (phenotype %in% names(ph.no.cov)) {
    print(aggregate(sfon$pheno[phenotype], by=list(sfon$geno[[chr]]$data[,mname]), FUN=mean, na.rm=T))
  }
}



#Confirming Actual Genotypes New#
sfon$geno[[4]]$data[,"7187"]
sfon$pheno$Fish.ID


###3J Specific chromosome-wide GxS interaction effect queries ####
# Save out the Additional_file_2 as a .csv to find out which are significant
add.file.2 <- read.csv(file = "02_data/additional_file_2_2017-01-05.csv")

# for the current version, need columns
names(add.file.2)
keep <- c("phenotype", "sex.cov.", "chr", "mname")
add.file.2 <- add.file.2[, keep]
head(add.file.2)
add.file.2 <- add.file.2[complete.cases(add.file.2),] # remove any rows with NA (this is to get rid of sex cov rows)
head(add.file.2)
length(add.file.2$mname)

# what phenos are to be used
unique(x = add.file.2$phenotype)
add.file.2$phenotype

# these are the actual names
names(ph.yes.cov)
names(ph.no.cov)
ph.sex.sp

# phenos
phenos.w.cov <- c(  rep("weight.g_0509", times = 3)
                  , rep("weight.g_0709", times = 3)
                  , rep("weight.g_1109", times = 2)
                  , rep("leng.cm_0509", times = 4)
                  , rep("leng.cm_0709", times = 5)
                  , rep("leng.cm_1109", times = 3)
                  , "hematocr"
                  , "cort.delta"
                  , "chlor.delta"
                  , "osmo.delta"
                  )
phenos.w.cov

chr.cols.w.cov <- c(1:20,30:34)
chr.w.cov <- add.file.2$chr[chr.cols.w.cov]

# double check
add.file.2$phenotype[chr.cols.w.cov]
phenos.w.cov # should match, but need to use this one b/c of changed pheno names


# Todo#
length(chr.cols.w.cov) # why is this wrong?


# CORTISOL WITH INTERACTION EFFECT
# test with a single pheno/chr combination
out.am.cort.delta.6 <- scanone(sfon, method = "hk", addcov = sex,
                             pheno.col = "cort.delta", chr = 6)
out.im.cort.delta.6 <- scanone(sfon, method = "hk", addcov = sex,
                             pheno.col = "cort.delta", chr = 6
                             , intcovar=sex)
# run permutations
num.perm.chr <- 10000
num.cluster <- 3
set.seed(54955149)
operm.am.cort.delta.6 <- scanone(sfon, addcovar=sex, n.perm=num.perm.chr,
                    pheno.col = "cort.delta", chr = 6,
                    n.cluster = num.cluster)
set.seed(54955149)
operm.im.cort.delta.6 <- scanone(sfon, addcovar=sex, n.perm=num.perm.chr,
                                 pheno.col = "cort.delta", chr = 6,
                                 n.cluster = num.cluster
                                 , intcovar=sex)
# analyze
head(operm.am.cort.delta.6) #Note: no pvals assoc. w scanone obj, pvals intro'd by 'summary'
out.im.cort.delta.6[1,] # each model obj. contains LOD scores, which can be compared
out.am.cort.delta.6[1,]
out.im.cort.delta.6[1,] - out.am.cort.delta.6[1,]
## for easy use, combine full and interactive models
out.fmim.cort.delta.6 <- c(out.im.cort.delta.6, out.im.cort.delta.6 - out.am.cort.delta.6, 
              labels=c("f","i"))
names(out.fmim.cort.delta.6) # note the positions of the full and interactive values

## for easy use, combine full and interactive models with the permutations applied
operm.fmim.cort.delta.6 <- cbind(operm.im.cort.delta.6, operm.im.cort.delta.6 - operm.am.cort.delta.6,
                    labels=c("f","i"))
colnames(operm.im.cort.delta.6) # NOTE that this contains all phenotypes, including those that do not require sex as cov.
colnames(operm.fmim.cort.delta.6)

# COMPARE RESULTS TO IDENTIFY SIG INTERACTION #
summary(operm.fmim.cort.delta.6, alpha=0.05) # provides the LOD threshold for p ≤ 0.05

# Check for significance in the full model and in the interactive model alone
summary(out.fmim.cort.delta.6, perms=operm.fmim.cort.delta.6, format="tabByCol"
        #, alpha=0.15
        , pvalues=TRUE
        )
# If high LOD comes from the interaction, rather than additive portion, sig ifx.
# If trait has no sig ifx, check sig in additive model only (interaction uses up power)
summary(out.am.cort.delta.6, perms=operm.am.cort.delta.6, format="tabByCol", alpha=0.1, pvalues=T)

# in conclusion, 
# no significant interaction effect
# lod.f p = 0.0364
# lod.i p = 0.177



# Osmolality WITH INTERACTION EFFECT on chr 40
# test with a single pheno/chr combination
out.am.osmo.delta.40 <- scanone(sfon, method = "hk", addcov = sex,
                               pheno.col = "osmo.delta", chr = 40)
out.im.osmo.delta.40 <- scanone(sfon, method = "hk", addcov = sex,
                               pheno.col = "osmo.delta", chr = 40
                               , intcovar=sex)
# run permutations
set.seed(54955149)
operm.am.osmo.delta.40 <- scanone(sfon, addcovar=sex, n.perm=num.perm.chr,
                                 pheno.col = "osmo.delta", chr = 40,
                                 n.cluster = num.cluster)
set.seed(54955149)
operm.im.osmo.delta.40 <- scanone(sfon, addcovar=sex, n.perm=num.perm.chr,
                                 pheno.col = "osmo.delta", chr = 40,
                                 n.cluster = num.cluster
                                 , intcovar=sex)
# analyze
head(operm.am.osmo.delta.40) #Note: no pvals assoc. w scanone obj, pvals intro'd by 'summary'
out.im.osmo.delta.40[1,] # each model obj. contains LOD scores, which can be compared
out.am.osmo.delta.40[1,]
out.im.osmo.delta.40[1,] - out.am.osmo.delta.40[1,]
## for easy use, combine full and interactive models
out.fmim.osmo.delta.40 <- c(out.im.osmo.delta.40, out.im.osmo.delta.40 - out.am.osmo.delta.40, 
                           labels=c("f","i"))
names(out.fmim.osmo.delta.40) # note the positions of the full and interactive values

## for easy use, combine full and interactive models with the permutations applied
operm.fmim.osmo.delta.40 <- cbind(operm.im.osmo.delta.40, operm.im.osmo.delta.40 - operm.am.osmo.delta.40,
                                 labels=c("f","i"))
colnames(operm.im.osmo.delta.40) # NOTE that this contains all phenotypes, including those that do not require sex as cov.
colnames(operm.fmim.osmo.delta.40)

# COMPARE RESULTS TO IDENTIFY SIG INTERACTION #
summary(operm.fmim.osmo.delta.40, alpha=0.05) # provides the LOD threshold for p ≤ 0.05

# Check for significance in the full model and in the interactive model alone
summary(out.fmim.osmo.delta.40, perms=operm.fmim.osmo.delta.40, format="tabByCol"
        #, alpha=0.15
        , pvalues=TRUE
)
# If high LOD comes from the interaction, rather than additive portion, sig ifx.
# If trait has no sig ifx, check sig in additive model only (interaction uses up power)
summary(out.am.osmo.delta.40, perms=operm.am.osmo.delta.40, format="tabByCol", alpha=0.1, pvalues=T)

# in conclusion, osmo 
# no significant interaction effect 
# lod.f = 0.027
# lod.i = 0.365


# weight.g_0709 WITH INTERACTION EFFECT (chr20)
# test with a single pheno/chr combination
out.am.weight.g_0709.20 <- scanone(sfon, method = "hk", addcov = sex,
                                pheno.col = "weight.g_0709", chr = 20)
out.im.weight.g_0709.20 <- scanone(sfon, method = "hk", addcov = sex,
                                pheno.col = "weight.g_0709", chr = 20
                                , intcovar=sex)
# run permutations
set.seed(54955149)
operm.am.weight.g_0709.20 <- scanone(sfon, addcovar=sex, n.perm=num.perm.chr,
                                  pheno.col = "weight.g_0709", chr = 20,
                                  n.cluster = num.cluster)
set.seed(54955149)
operm.im.weight.g_0709.20 <- scanone(sfon, addcovar=sex, n.perm=num.perm.chr,
                                  pheno.col = "weight.g_0709", chr = 20,
                                  n.cluster = num.cluster
                                  , intcovar=sex)
# analyze
head(operm.am.weight.g_0709.20) #Note: no pvals assoc. w scanone obj, pvals intro'd by 'summary'
out.im.weight.g_0709.20[1,] # each model obj. contains LOD scores, which can be compared
out.am.weight.g_0709.20[1,]
out.im.weight.g_0709.20[1,] - out.am.weight.g_0709.20[1,]
## for easy use, combine full and interactive models
out.fmim.weight.g_0709.20 <- c(out.im.weight.g_0709.20, out.im.weight.g_0709.20 - out.am.weight.g_0709.20, 
                            labels=c("f","i"))
names(out.fmim.weight.g_0709.20) # note the positions of the full and interactive values

## for easy use, combine full and interactive models with the permutations applied
operm.fmim.weight.g_0709.20 <- cbind(operm.im.weight.g_0709.20, operm.im.weight.g_0709.20 - operm.am.weight.g_0709.20,
                                  labels=c("f","i"))
colnames(operm.im.weight.g_0709.20) # NOTE that this contains all phenotypes, including those that do not require sex as cov.
colnames(operm.fmim.weight.g_0709.20)

# COMPARE RESULTS TO IDENTIFY SIG INTERACTION #
summary(operm.fmim.weight.g_0709.20, alpha=0.05) # provides the LOD threshold for p ≤ 0.05

# Check for significance in the full model and in the interactive model alone
summary(out.fmim.weight.g_0709.20, perms=operm.fmim.weight.g_0709.20, format="tabByCol"
        #, alpha=0.15
        , pvalues=TRUE
)

# no significant interaction effect here either.
# lod.f pval = 0.042
# lod.i pval = 0.416


####### IDENTIFY ACTUAL GENOS - NEEDS WORK!! ######
# this should all be correct but required manual back-calc from Rqtl to joinmap to haplotypes from STACKS
# and also manual obtaining of the column for the required marker

#check out a semi-informative marker
par(mfrow=c(1,2))
marker.count <- which(names(sfon$geno[[17]]$data[1,]) == "19336") #find a marker by name on chr 17
effectplot(sfon, mname1="19336",  mname2="sex",
           mark1=sfon$geno[[17]]$data[,marker.count], 
           geno1=c("CC", "CT"), # may or may not be correct
           mark2=sex, geno2=c("F","M"), 
           pheno.col = "leng.cm_0709",  #calling by name is safer
           #add.legend=F, 
           ylim = c(18,30),
           ylab = "length (cm) at T2",
           main = "",
           xlab = "",
           #add.legend = F
) 
#legend("topleft", legend = c("CC","CT"), fill = c("red","blue"))
# may want to add legend separately when generating figs for paper to choose position

# extra investigation, confirm that males larger than females
boxplot(pheno.df$leng.cm_0709 ~ pheno.df$sex, las = 1, 
        ylim = c(18,30)
        )  

# extra investigation, confirm direction of genotype x phenotype interaction
rqtl.19336 <- sfon$geno[[17]]$data[,143]
rqtl.19336.pheno <- sfon$pheno$leng.cm_0709
boxplot(rqtl.19336.pheno ~ rqtl.19336, las = 1)


# semi-informative marker on chr 15:
#example with covariate, but using the actual data, not imputed data
marker.count <- which(names(sfon$geno[[15]]$data[1,]) == "66075") #at position 100
effectplot(sfon, mname1="66075",  mname2="Sex",
           mark2=sex, geno2=c("F","M"), pheno.col = "leng.cm_0709", 
           mark1=sfon$geno[[15]]$data[,marker.count],
           add.legend=T, ylim = c(18, 30), main = ""
)
#legend("topleft", legend = c("AC","BC"), fill = c("black","red"))


######### Appendix 1. FORMULA FOR PLOTTING LOD CURVES #######
par(mfrow=c(1,1), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

#find your LOD column (start counting after chr and pos)
names(out.am) # needs covariate
names(all.out.0.nocov) # not needing covariate

# first designate scanone model to plot
# choose from 
#no covar model: all.out.0.nocov gwperms: all.out.0.nocov.perm
#additive covar: out.am gwperms: operm.am
#interactive: model: out.fmim; perms: operm.fmim NOTE: contains both the full and isolated interactive

# start by setting your variables
POI <- "weight.g_1109" #set your pheno
scanone.mod <- out.am
scanone.perms <- operm.am

lod.col <- which(names(scanone.mod) == POI) - 2 #-2 is to account for chr and pos
lod.col2 <- which(names(scanone.mod) == POI2) - 2 #-2 is to account for chr and pos


# then plot
plot(scanone.mod, alternate.chrid=T, lodcolumn=lod.col,
     #ylim = c(0,8), 
     ylab = POI,
     xlab = "",
     bandcol="gray70")
abline(h=summary(scanone.perms[,POI], 0.05), lty=1)


######### Appendix 2. PLOT QTL ON GENETIC MAP #######
qtl.pos <- c(28.3,198,162,83,115,261,169,112,89.9,46.9,150,35.7,139,31.4,66.8,47.8,114,43.3,60.3,138,22.6,139) #positions for 
qtl.chr <- c(4,5,20,3,4,5,20,34,16,39,20,3,9,36,5,24,6,38,40,24,4,25)
qtl.name <- c(rep("weight", times = 3), rep("length", times=5), "cfact1", "cfact2", "cfact3",
              rep("sgr2-3", times=3), "egg.diam", "sperm.diam", "d.cort", "d.chlor", "d.osmo", "ghr",
              rep("hemato", times=2))
#qtl.name <- paste("QTL", seq(1:length(qtl.pos)), sep = "")
#qtl.name <- seq(1:length(qtl.pos))


all.qtl <- makeqtl(cross=sfon, chr=qtl.chr, pos = qtl.pos, qtl.name=qtl.name, what = "draws")

par(mfrow=c(1,1), mar= c(3.5,3.5,1,1.5) + 0.2, mgp = c(2.5,0.75,0))
plot(all.qtl, chr = qtl.chr, alternate.chrid = T
     , col = c(rep("blue", times=14), rep("black", times=2), rep("red", times=6))
     #, horizontal = T
     #, justdots=T
      , main = ""
)

# save as 8x5
