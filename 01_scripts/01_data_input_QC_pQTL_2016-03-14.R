# Analysis of the phenotypic QTL from Sauvage et al with the new map from Sutherland et al.
# analysis v2.0
# Part 01/03 (data input and quality control) in R/QTL (cp cross type)

# Clean the space
rm(list=ls())

# Load R/qtl
#install.packages("qtl")
require(qtl)

# Set working directory
# setwd("~/Documents/bernatchez/01_Sfon_projects/03_Sfon_pQTL/sfon_pqtl")

##### IMPORT DATA #####
sfon <- read.cross(format="mapqtl"
                   , dir="02_data/", 
                   genfile = "Sfon_female_map_v4.3.loc", 
                   mapfile = "Sfon_female_map_v4.3.map", 
                   phefile = "SFQTL_phenotypes-full.qua",
                   genotypes = NULL, 
                   na.strings=c("NA","--"))
sfon <- jittermap(sfon)
orig.map <- pull.map(sfon)

plot.map(sfon, alternate.chrid = T, main ="", xlab = "Linkage group")
text(x = 35, y = 300, 
     labels = paste("n =", sum(nmar(sfon)), "markers") )


#####1A PHENOTYPE OUTLIER or TRANSFORMATION#####
nphe(sfon) #40 phenos
names(sfon$pheno)
names(sfon$pheno)[18] <- "hep.som.ind" # this was named incorrectly in raw data (need better way to do this)

#remove phenotypes not of interest
# NOTE: as an improvement, I think we should keep pre-stress values too for testing.
drop.phenos <- c(5,20,21,23,24,26,27) ## to fix, change this to the actual names that you want to remove (safer)
names(sfon$pheno[drop.phenos]) #display names of phenotypes to remove
sfon$pheno <- sfon$pheno[-drop.phenos]
names(sfon$pheno)
nphe(sfon) #33 retained (one is Fish.ID)

#plot raw phenotypes
par(mfrow=c(6,6), mar= c(2,3,1,1) + 0.2, mgp = c(2,0.75,0))
for(i in 2:33) 
  plot.pheno(sfon, pheno.col=i)
names(sfon$pheno)[29:33] #gene expr phenos are not normally distrib. they are prob linear.

# note that there are some fish with diff maturity index
par(mfrow=c(1,1))
plot.pheno(sfon, pheno.col="matur")
table(sfon$pheno$matur=="+++") # 2
table(sfon$pheno$matur=="++") # 15
table(sfon$pheno$matur=="-") # 8
boxplot(sfon$pheno$leng.cm_0709 ~ sfon$pheno$matur * sfon$pheno$sex, las = 1, ylab = "leng.cm_0709")
# and this has some indication for outliers in length.


# ASIDE # Keep in mind that there is variation in the data in terms of matur index
# May want to try full analysis with only samples with matur of "+"
sfon.matur <- subset(sfon, ind = which(sfon$pheno$matur=="+"))
plot.pheno(sfon.matur, pheno.col = "matur")
# Note the following
sex.mature <- as.numeric(pull.pheno(sfon.matur, "sex") == "M")
sex <- as.numeric(pull.pheno(sfon, "sex") == "M")
par(mfrow=c(1,2))

# note no effect for BC04 and BC05 in terms of broad peaks, but BC03 no longer elevated
plot(scanone(sfon, method="hk", addcov=sex, pheno.col="leng.cm_0709"), main = "sfon", ylab = "leng.cm_0709", alternate.chrid = T, ylim = c(0,6))
plot(scanone(sfon.matur, method="hk", addcov=sex.mature, pheno.col="leng.cm_0709"), main = "sfon.matur", ylab = "leng.cm_0709", alternate.chrid = T, ylim = c(0,6))

# note the significance viewed at BC38 is gone.
plot(scanone(sfon, method="hk", addcov=sex, pheno.col="chlor.delta"), main = "sfon", ylab = "chlor.delta", alternate.chrid = T, ylim = c(0,5))
plot(scanone(sfon.matur, method="hk", addcov=sex.mature, pheno.col="chlor.delta"), main = "sfon.matur", ylab = "chlor.delta", alternate.chrid = T, ylim = c(0,5))
# END ASIDE #

#log transform linear gene.expression
for(i in 29:33) {
  sfon$pheno[i] = log2(sfon$pheno[i])
}
#re-plot and can see more normal
nphe(sfon)
par(mfrow=c(6,6), mar= c(2,3,1,1) + 0.2, mgp = c(2,0.75,0))
for(i in 2:33) 
  plot.pheno(sfon, pheno.col=i)

## Remove outliers
# Outliers are visible for the following phenotypes:
outlier.containing.phenos <- c("TCS_T1.T2", "TCS_T2.T3", "leng.cm_0709"
                    , "condit.fact_T2", "osmo.delta", "male.sperm.diam"
                    )

# Loop across and remove outliers for each of the above phenotypes (3 x sd)
low.value = NULL; high.value = NULL; pheno.this.loop = NULL
for(pheno in 1:length(outlier.containing.phenos)) {
  pheno.this.loop <- outlier.containing.phenos[pheno]
  # find the > 3 * SD outlier range for the phenotype
  low.value <- mean(sfon$pheno[[pheno.this.loop]], na.rm = T) - 3*sd(sfon$pheno[[pheno.this.loop]], na.rm = T)
  high.value <- mean(sfon$pheno[[pheno.this.loop]], na.rm = T) + 3*sd(sfon$pheno[[pheno.this.loop]], na.rm = T)
  # change the outliers to NA
  sfon$pheno[[pheno.this.loop]][sfon$pheno[[pheno.this.loop]] < low.value 
                       | sfon$pheno[[pheno.this.loop]] > high.value] <- NA
  }


## Rename individuals with sex mistakenly identified (based on gonad measure and transcriptome)
sfon$pheno$sex[sfon$pheno$male.sperm.conc != "NA"] <- "M"
sfon$pheno$sex[sfon$pheno$fem.egg.diam != "NA"] <- "F"


#####1B SEGREGATION DISTORTION####
gt = geno.table(sfon) #check for segregation distortion
names(gt)
ld <- gt[gt$P.value < 0.01,] #identify markers with seg.dis p < 0.01
length(ld$P.value) #157 markers with segregation distortion
summary(ld$chr) 
ld.count.chr <- summary.factor(ld$chr)
ld.count.chr[ld.count.chr > 13]  #chromosomes with many segregation distortion markers

sfon <- drop.markers(sfon, rownames(ld)) # remove the segregation distortion markers
gt2 = geno.table(sfon) #recheck for segregation distortion
ld2 <- gt2[gt2$P.value < 0.01,] 
ld2 #all have been removed


#####1C CHECK FOR IDENTICAL GENOTYPES#####
cg = comparegeno(sfon)
hist(cg, breaks=200, xlab="Proportion of identical genotypes", main = "Sfon map",
     xlim = c(0:1)) #result: no outliers in curve of identical genotypes


#####1D ESTIMATE RECOMBINATION FRACTION########
sfon <- est.rf(sfon)
#plot.rf(sfon, what = "both", main = "Sfon (RF/LOD)", alternate.chrid=T) #long plot


#####1E COUNT CROSSOVER EVENTS #####
nxo <- countXO(sfon) #check for number of obligate crossovers by individual
plot(nxo, ylab="No. crossovers", xlab = "order")  
sfon$pheno$Fish.ID[nxo>300] #The 150th individual has 1093 crossovers, whereas everyone else has ~50-200
text("SFO-171", x= 100, y=1050)

sfon <- subset(sfon, ind = -150) #remove individual with abnormally high crossovers

# Can also take a look at how many crossovers are needed per chromosome
colSums(countXO(sfon, bychr = T))


#####1F MARKER ORDER AND MAP PLOT #####
# to be done on final map only
#check for correct ordering of markers:
#rip <- ripple(sfon, chr=33) #xxx total orders
#summary(rip) #gives the initial (input) marker order; then the order w lowest number of obligate crossovers
#LG01.switch <- switch.order(LG01.imp, chr="LG01.imputed", order=rip[2,]) #to properly do, need the chromosome
#nm <- est.map(sfon, error.prob=0.001) ##estimate map using existing genodaata

#plot map
par(mfrow=c(1,1))
plot.map(sfon, alternate.chrid=T, main = "")

#####1G MISSING GENOTYPES ######
plot.info(sfon, alternate.chrid = T, fourwaycross = "all") #check for missing data
plot.info(sfon, alternate.chrid = T, fourwaycross = "AB") # missing alleles from parent 1
plot.info(sfon, alternate.chrid = T, fourwaycross = "CD") # missing alleles from parent 2
# note that a larger amount of missing information for parent 1 (the male)
# this only occurs due to the male-specific markers not being included

z <- plot.info(sfon, alternate.chrid = T, step=0) #step gives density of grid; #step = 0; calc only performed at the markers
dim(z) #5997 rows, three cols
head(z) #name.marker; chr; pos; misinfo.entropy
#z[z[,1]==14,] #looks at chr14 specifically

# missing genotypes per marker
par(mar = c(3.4,3.4,2,1), mfrow = c(1,2))
hist(nmissing(sfon, what="mar"), 
     breaks = 50, main = "", las = 1,
     xlab = "Number of missing genotypes per marker") #counts # of missing genos for each marker
abline(v = 5)
hist(nmissing(sfon, what="ind"), 
     breaks = 50, main = "", las = 1,
     xlab = "Number of missing genotypes per indiv") #counts # of missing genos for each marker
abline(v = 100)

# Bit of exploration of which samples/markers are missing high numbers of genotypes
marker.miss.geno <- nmissing(sfon, what="mar") 
marker.miss.geno[marker.miss.geno>75]

indiv.miss.geno <- nmissing(sfon, what="ind")
sfon$pheno$Fish.ID[which(indiv.miss.geno > 2000)]

## to do outside of R to find out what types of markers have highest missing data (most are nnxnp, but most are..)
# write.csv(test, file = "03_results/markers_missing_more_than_75_indiv_genos.csv", quote = F)
# awk -F, '{ print $1 }' 03_results/markers_missing_more_than_75_indiv_genos.csv | grep -vE '$^' > 03_results/markers_missing_more_than_75_indiv_genos_clean.txt
# grep -f 03_results/markers_missing_more_than_75_indiv_genos_clean.txt 02_data/*.loc | less

#####1H GENOTYPE ERRORS #####
# Takes a very long time, use multiple cores, and only do with the final version of the map
#newmap <- est.map(sfon, error.prob=0.01) #estimate map
#hyper <- replace.map(hyper, newmap) #replacing map with the estimated one from data
#test <- calc.errorlod(sfon)
#top <- top.errorlod(test, cutoff=5) #above a specified cutoff
#plotErrorlod
#top #note, ID column here is numeric index. if pheno named 'id' it would have been used

#plot.geno(hyper, 16, top$id[top$chr==16], cutoff=5)
#a small number of genotyping errors will not have much an influence
# on results - but if you see a QTL in the region, may want to take a second
# look at raw data.


#####1I HMM CALCULATE GENOPROBABILITES####
sfon <- calc.genoprob(sfon, step = 0, error.prob =0.0001) #hmm to calc prob of true underlying genos given data


#####1J WHICH PHENOs REQ SEX AS COVARIATE?#####
names(sfon$pheno)
pheno.df <- data.frame(sfon$pheno[,-1]) #drop Fish.ID
row.names(pheno.df) <- sfon$pheno[,1] #use Fish.ID as row.names
names(pheno.df)

# collect the phenos to test for an association with sex
phenos <- c("weight.g_0509", "weight.g_0709", "weight.g_1109",
            "TCS_T1.T2", "TCS_T2.T3", "TCS_T1.T3", "leng.cm_0509","leng.cm_0709","leng.cm_1109",
            "condit.fact_T1", "condit.fact_T2", "condit.fact_T3", "weight_liver.g", "hep.som.ind", 
            "hematocr", "cort.delta", "osmo.delta", "chlor.delta",
            "plasma.chlor", "plasma.osmo", "plasma.gluc", "hep.glyco",
            "ghr", "igf1", "igfr1", "ef1","bactin") #collect phenos of interest; do not include sex-specific phenos
length(phenos) #27 indiv phenos of interest to compare between sexes
names(sfon$pheno[, phenos]) 

# boxplot phenotypes including covariate
par(mfrow=c(6,5),  mar= c(2,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
for(i in phenos)
  boxplot(pheno.df[,i] ~ pheno.df$sex, main = names(pheno.df[i]))

# test for signif effect of covariate
pvals.pxs <- NULL
for(i in phenos) {
  pvals.pxs = c(pvals.pxs, anova(aov(pheno.df[ ,i] ~ pheno.df$sex))$`Pr(>F)`[1])
}
names(pvals.pxs) <- phenos
pvals.pxs[pvals.pxs < 0.2] # this displays which phenotypes could use a covariate
length(pvals.pxs[pvals.pxs < 0.2]) # number of phenos needing covar
ph.yes.cov <- which(pvals.pxs < 0.2) # phenos needing covariate
#ph.yes.cov.num <- match(names(ph.yes.cov), names(pheno.df))

ph.no.cov <- which(pvals.pxs >= 0.2) # phenos NOT needing covariate
#ph.no.cov.num <- match(names(ph.no.cov), names(pheno.df))

ph.sex.sp <- c("fem.egg.diam","male.sperm.conc","male.sperm.diam")


#####1K CORRELATION OF PHENOTYPES #####
library('corrplot') #package corrplot
names(pheno.df)
to.cor.phenos <- c("weight.g_0709",  "leng.cm_0709", "TCS_T1.T3"
                   , "condit.fact_T2", "hep.som.ind", "weight_liver.g"
                   , "chlor.delta", "cort.delta", "osmo.delta"
                   , "fem.egg.diam", "male.sperm.conc", "male.sperm.diam"
                   , "hematocr", "plasma.chlor", "plasma.osmo", "plasma.gluc"
                   , "hep.glyco", "ghr", "igf1", "igfr1")
# could also do that with subtracting the ones not wanted

cor.set <- cor(pheno.df[,to.cor.phenos]
    , use = "pairwise.complete.obs" # cor bw ea pair var is comp using all compl pairs of obs on those var # note: only works with "pearson" method    
    )

par(mfrow=c(1,1), mar= c(3,5,3,1) + 0.2, mgp = c(2,0.75,0))
corrplot(cor.set, method = "circle"
         #, order="hclust"
         , addshade = "all"
         #, type = "lower"
         #, outline = T
         , tl.col = "black"
         , tl.cex = 0.9) #plot matrix

#problem with missing data where all are missing between phenos.
?corrplot



#### NEEDS TO BE FIXED ####
######### Find phenotype averages and standard deviations ####
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


#######NOW GO TO PART 2 SCRIPT TO DO SCANONE PERM TESTS#######
#export CLEANED INPUT DATA FOR SCANONE PERM TESTS #######
save.image(file = "02_data/sfon_01_output.RData") # save out existing data so that it can be reloaded without running all again

