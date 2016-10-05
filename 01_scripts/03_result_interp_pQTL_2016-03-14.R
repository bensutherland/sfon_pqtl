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
#load("02_data/sfon_01_output.RData")
#load("02_data/sfon_02_output.RData")
#load("02_data/sfon_02_output_gw_and_chr.RData") # both chromosome and genome-wide significance
#load("02_data/sfon_02_output_chromosome-wide_p0.01_1000perms.RData") # 
#load("02_data/sfon_02_output_chromosome-wide_p0.01_10000perms.RData") # 

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


######## EXTRACT RESULTS OF CHR-WIDE SIG ######
summary(scanone.mods[["weight.g_0509_1"]], threshold = 1.2) # for example
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

# Try again:
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

write.csv(chr.wide.output, file = "chr_wide_p0.01_w_10000perms.csv", quote = F, row.names=F)
# then in terminal:
# grep -A1 'lod' test.csv | grep -v '--' - | awk -F, '{ print $2 }' | less

# if want it back in:
test <- read.csv(file = "test.csv", header = F)




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


########## PVE Genome Wide ######
# for this use the following formulae (per trait)
?makeqtl()
qtl.trait <- makeqtl(sfon, chr=c(), pos=c())
?fitqtl()
out.fq <- fitqtl(sfon, qtl=qtl.trait, formula=y~Q1+Q2+Q3) # put all of the QTL in one formula
summary(out.fq) # this will give the PVE
plot(qtl.object) # this will put the QTL on your genetic map

# example with weight.g_0509
sex.df <- as.data.frame(sex)
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
qtl.weight.g_0709_sex <- makeqtl(sfon, chr=c("4","5"), pos = c(28.3, 198.5), what = "draws")
out.weight.g_0709_sex <- fitqtl(sfon, pheno.col="weight.g_0709", qtl=qtl.weight.g_0709_sex
                                , formula=y~Q1+Q2, covar=sex)
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

##### PVE chromosome-wide #######
# TCS_T1.T2
qtl.chr.TCS_T1.T2 <- makeqtl(sfon, chr=c("7","9","15","27","40"), pos = c(59.1, 184, 77.7, 121, 95.6), what = "draws")
out.chr.TCS_T1.T2 <- fitqtl(sfon, pheno.col="TCS_T1.T2", qtl=qtl.chr.TCS_T1.T2, formula=y~Q1+Q2+Q3+Q4+Q5)
summary(out.chr.TCS_T1.T2)

# TCS_T2.T3
qtl.chr.TCS_T2.T3 <- makeqtl(sfon, chr=c("3","9","23","36","37"), pos = c(35.7, 139, 131, 31.4, 24.3), what = "draws")
out.chr.TCS_T2.T3 <- fitqtl(sfon, pheno.col="TCS_T2.T3", qtl=qtl.chr.TCS_T2.T3, formula=y~Q1+Q2+Q3+Q4+Q5)
summary(out.chr.TCS_T2.T3)

# TCS_T1.T3
qtl.chr.TCS_T1.T3 <- makeqtl(sfon, chr=c("3","15","36"), pos = c(35.7, 170, 61.8), what = "draws")
out.chr.TCS_T1.T3 <- fitqtl(sfon, pheno.col="TCS_T1.T3", qtl=qtl.chr.TCS_T1.T3, formula=y~Q1+Q2+Q3)
summary(out.chr.TCS_T1.T3)

# condit.fact_T1
qtl.chr.condit.fact_T1 <- makeqtl(sfon, chr=c("16","24","36"), pos = c(89.9, 79.8, 88.4), what = "draws")
out.chr.condit.fact_T1 <- fitqtl(sfon, pheno.col="condit.fact_T1", qtl=qtl.chr.condit.fact_T1, formula=y~Q1+Q2+Q3)
summary(out.chr.condit.fact_T1)

# condit.fact_T2
qtl.chr.condit.fact_T2 <- makeqtl(sfon, chr=c("1","3","14","19","20","24","39"), 
                                  pos = c(266,6.78,97.9,122,51.1,67.9,46.9), what = "draws")
out.chr.condit.fact_T2 <- fitqtl(sfon, pheno.col="condit.fact_T2", qtl=qtl.chr.condit.fact_T2, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7)
summary(out.chr.condit.fact_T2)

# condit.fact_T3
qtl.chr.condit.fact_T3 <- makeqtl(sfon, chr=c("6","19","20","24","38"), 
                                  pos = c(19.7,7.37,150,79.8,46.4), what = "draws")
out.chr.condit.fact_T3 <- fitqtl(sfon, pheno.col="condit.fact_T3", qtl=qtl.chr.condit.fact_T3, formula=y~Q1+Q2+Q3+Q4+Q5)
summary(out.chr.condit.fact_T3)

# chlor.delta
qtl.chr.chlor.delta <- makeqtl(sfon, chr=c("38"), pos = c(43.3), what = "draws")
out.chr.chlor.delta <- fitqtl(sfon, pheno.col="chlor.delta", qtl=qtl.chr.chlor.delta, formula=y~Q1)
summary(out.chr.chlor.delta)

# plasma.chlor
qtl.chr.plasma.chlor <- makeqtl(sfon, chr=c("3","9","20","37"), pos = c(174,200,175,24.3), what = "draws")
out.chr.plasma.chlor <- fitqtl(sfon, pheno.col="plasma.chlor", qtl=qtl.chr.plasma.chlor, formula=y~Q1+Q2+Q3+Q4)
summary(out.chr.plasma.chlor)

# plasma.gluc
qtl.chr.plasma.gluc <- makeqtl(sfon, chr=c("17","36"), pos = c(82.5,31.4), what = "draws")
out.chr.plasma.gluc <- fitqtl(sfon, pheno.col="plasma.gluc", qtl=qtl.chr.plasma.gluc, formula=y~Q1+Q2)
summary(out.chr.plasma.gluc)

# fem.egg.diam
qtl.chr.fem.egg.diam <- makeqtl(sfon, chr=c("5","19"), pos = c(66.8,105), what = "draws")
out.chr.fem.egg.diam <- fitqtl(sfon, pheno.col="fem.egg.diam", qtl=qtl.chr.fem.egg.diam, formula=y~Q1+Q2)
summary(out.chr.fem.egg.diam)

# male.sperm.conc
qtl.chr.male.sperm.conc <- makeqtl(sfon, chr=c("13","22"), pos = c(25.9,38.8), what = "draws")
out.chr.male.sperm.conc <- fitqtl(sfon, pheno.col="male.sperm.conc", qtl=qtl.chr.male.sperm.conc, formula=y~Q1+Q2)
summary(out.chr.male.sperm.conc)

# male.sperm.diam
qtl.chr.male.sperm.diam <- makeqtl(sfon, chr=c("5","7","24","28","40"), 
                                  pos = c(243,221,47.8,59.5,20.9), what = "draws")
out.chr.male.sperm.diam <- fitqtl(sfon, pheno.col="male.sperm.diam", qtl=qtl.chr.male.sperm.diam, formula=y~Q1+Q2+Q3+Q4+Q5)
summary(out.chr.male.sperm.diam)


# weight.g_0509
qtl.chr.weight.g_0509 <- makeqtl(sfon, chr=c("4","5","9","20","27","34","36"), 
                                  pos = c(28.3,261,13.2,162,132,112,14.3), what = "draws")
out.chr.weight.g_0509 <- fitqtl(sfon, pheno.col="weight.g_0509", qtl=qtl.chr.weight.g_0509, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7)
summary(out.chr.weight.g_0509)

# weight.g_0709
qtl.chr.weight.g_0709 <- makeqtl(sfon, chr=c("4","5","20","24","27","36","39"), 
                                 pos = c(28.3,198,162,79.8,44.2,106,85.9), what = "draws")
out.chr.weight.g_0709 <- fitqtl(sfon, pheno.col="weight.g_0709", qtl=qtl.chr.weight.g_0709, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7)
summary(out.chr.weight.g_0709)

# weight.g_1109
qtl.chr.weight.g_1109 <- makeqtl(sfon, chr=c("3","5","9","20","24","28","34"), 
                                 pos = c(52.9,134,13.2,169,79.8,8.62,47.1), what = "draws")
out.chr.weight.g_1109 <- fitqtl(sfon, pheno.col="weight.g_1109", qtl=qtl.chr.weight.g_1109, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7)
summary(out.chr.weight.g_1109)


# leng.cm_0509
qtl.chr.leng.cm_0509 <- makeqtl(sfon, chr=c("5","10","12","16","19","20","27","30","34","39","40"), 
                                 pos = c(261,76.3,79.3,20.7,156,162,45.1,117,110,111,95.6), what = "draws")
out.chr.leng.cm_0509 <- fitqtl(sfon, pheno.col="leng.cm_0509", qtl=qtl.chr.leng.cm_0509, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q11)
summary(out.chr.leng.cm_0509)

# leng.cm_0709
qtl.chr.leng.cm_0709 <- makeqtl(sfon, chr=c("2","3","4","5","9","19","20","27","34","36","39"), 
                                 pos = c(260,83,115,261,13.2,134,169,45.1,112,14.3,0), what = "draws")
out.chr.leng.cm_0709 <- fitqtl(sfon, pheno.col="leng.cm_0709", qtl=qtl.chr.leng.cm_0709, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10+Q11)
summary(out.chr.leng.cm_0709)

# leng.cm_1109
qtl.chr.leng.cm_1109 <- makeqtl(sfon, chr=c("3","5","20","24","31","34","36","40"), 
                                 pos = c(52.9,134,169,79.8,126,47.1,93.5,0), what = "draws")
out.chr.leng.cm_1109 <- fitqtl(sfon, pheno.col="leng.cm_1109", qtl=qtl.chr.leng.cm_1109, formula=y~Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8)
summary(out.chr.leng.cm_1109)

# weight_liver.g
qtl.chr.weight_liver.g <- makeqtl(sfon, chr=c("1","2","35"), 
                                 pos = c(148,318,46.5), what = "draws")
out.chr.weight_liver.g <- fitqtl(sfon, pheno.col="weight_liver.g", qtl=qtl.chr.weight_liver.g, formula=y~Q1+Q2+Q3)
summary(out.chr.weight_liver.g)

# hep.som.ind
qtl.chr.hep.som.ind <- makeqtl(sfon, chr=c("3","35","36"), 
                                  pos = c(39.1,106,31.4), what = "draws")
out.chr.hep.som.ind <- fitqtl(sfon, pheno.col="hep.som.ind", qtl=qtl.chr.hep.som.ind, formula=y~Q1+Q2+Q3)
summary(out.chr.hep.som.ind)

# hematocr
qtl.chr.hematocr <- makeqtl(sfon, chr=c("3","4","25"), 
                               pos = c(237,22.6,139), what = "draws")
out.chr.hematocr <- fitqtl(sfon, pheno.col="hematocr", qtl=qtl.chr.hematocr, formula=y~Q1+Q2+Q3)
summary(out.chr.hematocr)

# cort.delta
qtl.chr.cort.delta <- makeqtl(sfon, chr=c("3","6","24"), 
                            pos = c(51.1,114,160), what = "draws")
out.chr.cort.delta <- fitqtl(sfon, pheno.col="cort.delta", qtl=qtl.chr.cort.delta, formula=y~Q1+Q2+Q3)
summary(out.chr.cort.delta)

# osmo.delta
qtl.chr.osmo.delta <- makeqtl(sfon, chr=c("5","40"), 
                              pos = c(9.2,60.3), what = "draws")
out.chr.osmo.delta <- fitqtl(sfon, pheno.col="osmo.delta", qtl=qtl.chr.osmo.delta, formula=y~Q1+Q2)
summary(out.chr.osmo.delta)

# plasma.osmo
qtl.chr.plasma.osmo <- makeqtl(sfon, chr=c("19","25","38"), 
                              pos = c(176,10.2,54.4), what = "draws")
out.chr.plasma.osmo <- fitqtl(sfon, pheno.col="plasma.osmo", qtl=qtl.chr.plasma.osmo, formula=y~Q1+Q2+Q3)
summary(out.chr.plasma.osmo)

# ghr
qtl.chr.ghr <- makeqtl(sfon, chr=c("11","24","36"), 
                               pos = c(53,138,43.9), what = "draws")
out.chr.ghr <- fitqtl(sfon, pheno.col="ghr", qtl=qtl.chr.ghr, formula=y~Q1+Q2+Q3)
summary(out.chr.ghr)

# igf1
qtl.chr.igf1 <- makeqtl(sfon, chr=c("4","8","9","13","40"), 
                       pos = c(227,42.9,13.2,173,7.8), what = "draws")
out.chr.igf1 <- fitqtl(sfon, pheno.col="igf1", qtl=qtl.chr.igf1, formula=y~Q1+Q2+Q3+Q4+Q5)
summary(out.chr.igf1)

# igfr1
qtl.chr.igfr1 <- makeqtl(sfon, chr=c("13","24","40"), 
                        pos = c(178,61.1,5.07), what = "draws")
out.chr.igfr1 <- fitqtl(sfon, pheno.col="igfr1", qtl=qtl.chr.igfr1, formula=y~Q1+Q2+Q3)
summary(out.chr.igfr1)



##### PHENO BY GENO GRAPHS ####

