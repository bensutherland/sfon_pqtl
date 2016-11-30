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
# Use marker list to remove markers
sfon_remove <- read.csv(file = "02_data/sfon_efxeg_and_hkxhk_markers.csv", header = F)
sfon_remove <- as.numeric(sfon_remove)
length(sfon_remove) # 1208 markers to remove
sfon_only_nnxnp <- drop.markers(sfon, sfon_remove) # drop markers
rm(sfon_remove)
# note that some may already be removed due to seg distortion filter

summary(sfon_only_nnxnp) #2509 markers
plot.map(sfon_only_nnxnp, alternate.chrid = T,
         main = "", xlab = "Linkage Group")
text(x = 35, y = 260, labels = paste("n =", sum(nmar(sfon_only_nnxnp)), "markers"))

##export subset data
# save.image(file = "02_data/sfon_01_output_subset_only_nnxnp.RData")

#### REMOVE SPECIFIC INDIV SEX ####
# Identify male and female individuals
ind.males = c(sfon_only_nnxnp$pheno$sex=="M")
ind.females = c(sfon_only_nnxnp$pheno$sex=="F")

sfon.males <- subset(sfon_only_nnxnp, ind=ind.males)
sfon.females <- subset(sfon_only_nnxnp, ind=ind.females)

#export subset data
save.image(file = "02_data/sfon_01_output_subset_only_nnxnp_and_mf_sep.RData")
