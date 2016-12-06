# Subset R/qtl object for marker type or for individual sex
# this assumes you have made marker lists for the marker types hkxhk and efxeg to remove (see readme.md)
# by doing
# grep -E 'efxeg|hkxhk' Sfon_female_map_v4.3.loc | awk 'BEGIN{ORS=","}1 { print $1 }' | sed 's/,$//g' > sfon_efxeg_and_hkxhk_markers.csv

# rm(list=ls())

# install libraries
require(qtl)

# set working directory
setwd("~/Documents/bernatchez/01_Sfon_projects/03_Sfon_pQTL/sfon_pqtl")

# Load Part 1 results:
load("02_data/sfon_01_output.RData")


#### REMOVE SPECIFIC MARKER TYPES ####
# Choose which removal set
#sfon_remove <- read.csv(file = "02_data/sfon_efxeg_and_hkxhk_markers.csv", header = F) # keep only nnxnp
#sfon_remove <- read.csv(file = "02_data/sfon_nnxnp_and_hkxhk_markers.csv", header = F) # keep only efxeg
#sfon_remove <- read.csv(file = "02_data/sfon_efxeg_markers.csv", header =F) # keep nnxnp and hkxhk
#sfon_remove <- read.csv(file = "02_data/sfon_nnxnp_markers.csv", header = F)
sfon_remove <- as.numeric(sfon_remove)

length(sfon_remove) # markers to remove
sfon_limited <- drop.markers(sfon, sfon_remove) # drop markers
rm(sfon_remove)
# note that some may already be removed due to seg distortion filter

summary(sfon_limited) # check markers
plot.map(sfon_limited, alternate.chrid = T,
         main = "", xlab = "Linkage Group")
text(x = 35, y = 260, labels = paste("n =", sum(nmar(sfon_limited)), "markers"))

# note that re-estimating recombination fraction after dropping those markers doesn't change map positions or total lengths etc.

# sfon <- est.rf(sfon) #doesn't seem to be necessary for crossover determination or for cM position

##export subset data
# save.image(file = "02_data/sfon_01_output_subset_only_nnxnp.RData")
# save.image(file = "02_data/sfon_01_output_subset_only_nnxnp_and_hkxhk.RData")
save.image(file = "02_data/sfon_01_output_subset_only_efxeg_and_hkxhk.RData") # to use for heterochiasmy



#### REMOVE SPECIFIC INDIV SEX ####
# Identify male and female individuals
ind.males = c(sfon_only_nnxnp$pheno$sex=="M")
ind.females = c(sfon_only_nnxnp$pheno$sex=="F")

sfon.males <- subset(sfon_only_nnxnp, ind=ind.males)
sfon.females <- subset(sfon_only_nnxnp, ind=ind.females)

#export subset data
save.image(file = "02_data/sfon_01_output_subset_only_nnxnp_and_mf_sep.RData")