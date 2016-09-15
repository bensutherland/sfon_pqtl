# Analysis of the phenotypic QTL from Sauvage et al with the new map from Sutherland et al.
# analysis v2.0
# Part 03/03 (result interpretation and figure generation)
# in R/QTL (cp cross type)

# rm(list=ls())

# install libraries
require(qtl)

# Set working directory
setwd("~/Documents/bernatchez/01_Sfon_projects/03_Sfon_pQTL/sfon_pqtl")

# Load data from permutation tests
load("02_data/sfon_01_output.RData")
load("02_data/sfon_02_output.RData")

# For most of the plotting (following) we need to simulate genotypes given observed marker data
sfon <- sim.geno(sfon, step=2.5,
                  error.prob=0.001, n.draws=256) #note that 4-way cross w sex-specific map assumes constant ratio of female:male recomb rates within inter-marker intervals
str((sfon$geno[[1]])) # note that there is now probabilities and draws added to each LG

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
     #ylim = c(0,8), 
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
summary(out.fmim, perms=operm.fmim, format="tabByCol", alpha=0.1, pvalues=TRUE)
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
abline(h=summary(scanone.perms[,POI], 0.05), lty=3)


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


########EXTRACT RESULTS OF CHR-WIDE SIG######
#e.g.
summary(scanone.mods[["weight.g_0509_1"]], threshold = 1.2)
summary(scanone.mods[["weight.g_0509_1"]], threshold = pheno.sig.lod.per.chr[1,"weight.g_0509"])

# Extract chromosome-wide significance from scanone object:
test <- NULL

test <- capture.output(
  for(pheno in all.phenos) {
    for(chr in selected.chrs) {
      print(c(pheno, chr), quote = F)
      print(summary(scanone.mods[[paste(pheno, "_", chr, sep="")]], threshold = pheno.sig.lod.per.chr[chr, pheno]), quote = F)
    }
  }
)
# can use gsub to remove the annoying quotes and backslash (NOT WORKING)
# gsub(pattern = "\"", replacement = "", x = test)
### this part will be in part 3 eventually ###


######### 3C RESULT INTERP, SEX AS A BINARY TRAIT #######
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

######### 3D PLOT DATA, sex markers #######
# find marker with highest LOD
find.marker(cross = sfon, chr = 35, pos = 109) 
# Note: if the marker is not fully informative, genotype plotting will vary each time due to imputation of genotypes
# so it may be more informative to find the closest fully informative marker as long as the lod is > thresh

plot.pxg(sfon, marker="43393", pheno.col = "binary") #infer=F doesn't prevent imputation         
# infer = T means that missing genos are filled with a random imputation and plotted in red

effectplot(sfon, pheno.col= "binary", mname1="83995") # still uses imputation


######### 3E PLOT DATA, EFFECT PLOTS #######

# This section needs some serious work #

# find a trait and marker with high lod, with an additive effect
summary(out.am, perms=operm.am, format="tabByCol", alpha=0.1, pvalues=T)
summary(out.fmim, perms=operm.fmim, format="tabByCol", alpha=0.1, pvalues=T) # make sure no interact
# e.g. for weight.g_0509 , consider chr 5 marker 85980 

par(mfrow=c(1,2), mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

#example with covariate, but using the actual data, not imputed data
# first lets find the actual data
sfon$geno[[5]]$data[1:5,1:5] # gives the genotype codes from the first five markers on chr 5 for the first five indiv (out of 170)
#Important## if using mark1 make sure to use geno1 or else arbitrary rqtl genotype labels will be used (i.e. first 2)

# as an example, check out a fully informative marker (with four distinct genotypes):
par(mfrow=c(1,1))
marker.count <- which(names(sfon$geno[[5]]$data[1,]) == "85980") #find a marker by name on chr 17
effectplot(sfon, mname1="128012",  
           mname2="Sex",
           mark1=sfon$geno[[5]]$data[,marker.count], 
           geno1=c("AG,AG","GA,AG","AG,AA","GA,AA"), #appropriate genos are added from stacks.markers.info.xlsx
           mark2=sex, geno2=c("F","M"), 
           pheno.col = "weight.g_0509", 
           #add.legend=F, 
           #ylim = c(20,28),
           ylab = "length (cm) at T2",
           main = "",
           xlab = ""
) 
#use geno1 to give the actual genotypes

# is there a significant chr-wide interaction


# automatic capture of the genotypes using 'mname' without assigning mark1 will
# use the genotype data that has been imputed, and thus may be slightly dependent on
# the sim.geno run.
effectplot(sfon, mname1 = "85980", mname2 = "Sex",
           mark2=sex, geno2 = c("F","M"),
           pheno.col = "weight.g_0509",
           add.legend=T,
           #ylim = c(20,28),
           main = "")


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


######### FORMULA FOR PLOTTING LOD CURVES #######
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


###### REMAINING QUESTIONS ########
###why is hep.som.ind so high for gw-sig ??
# note: so is TCS_T2.T3


## Exploring Phenotypes ##

######### preliminary phenotypes ####
names(sfon$pheno) # take 2:length(names(sfon$pheno))
length(names(sfon$pheno))
# find averages for the sexes
sex.sp.avg.df <- aggregate(sfon$pheno[,2:34], list(sfon$pheno$sex), na.rm = T, FUN=mean)
dim(sex.sp.avg.df)
sex.sp.avg.df

# and standard deviation
sex.sp.sd.df <- aggregate(sfon$pheno[,2:34], list(sfon$pheno$sex), na.rm = T, FUN=sd)
dim(sex.sp.sd.df)
sex.sp.sd.df

data <- rbind(sex.sp.avg.df, sex.sp.sd.df)
row.names(data) <- c("avg.F","avg.M","sd.F","sd.M") 

# because aggregate seems to not work with multiple values, will have to create a loop for that.
test <- NULL
# and length (sample sizes)
for(i in 2:34) {
  test <- aggregate(sfon$pheno[,i] ~ sfon$pheno$sex, FUN = length)
  print(names(sfon$pheno)[i])
  print(test)  
}

# To assign to object
test = NULL
test <- as.data.frame(c("temp", "temp"))
for(i in 2:34) {
  test <- cbind(test, aggregate(sfon$pheno[,i] ~ sfon$pheno$sex, FUN = length))
  #print(names(sfon$pheno)[i])
  #print(test)  
}
dim(test)

test[1:2,1:10]

sex.order <- test[1:2,2]

subset.samplesizes <- cbind(sex.order, 
                            test[, seq(3, ncol(test), by = 2)]
)
colnames(subset.samplesizes) <- c("Group.1", names(sfon$pheno)[-1])
rownames(subset.samplesizes) <- c("n.fem", "n.mal")
subset.samplesizes

# Bring together all data
dim(data)
names(data)
dim(subset.samplesizes)
names(subset.samplesizes)

avg.sd.n.df <- as.data.frame(rbind(data, subset.samplesizes))
str(avg.sd.n.df)
write.csv(x = t(avg.sd.n.df), file = "avg.sd.n.csv", quote = F)


########## PVE ######
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

# need to put in the covariate? doesn't seem to affect anything
summary(out.weight.g_0509.fq)

# example with weight.g_0709
qtl.weight.g_0709 <- makeqtl(sfon, chr=c("4","5"), pos = c(28.3, 198.5), what = "draws")
out.weight.g_0709 <- fitqtl(sfon, pheno.col="weight.g_0709", qtl=qtl.weight.g_0709, formula=y~Q1+Q2)
summary(out.weight.g_0709)

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

