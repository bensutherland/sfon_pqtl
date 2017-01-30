# Subset R/qtl object for marker type or for individual sex
# this assumes you have made marker lists for the marker types hkxhk and efxeg to remove (see readme.md)
# by doing
# 

# rm(list=ls())

# install libraries
require(qtl)

# set working directory
setwd("~/Documents/bernatchez/01_Sfon_projects/03_Sfon_pQTL/sfon_pqtl")

# Load Part 1 results:
load("02_data/sfon_01_output.RData")


#### REMOVE SPECIFIC MARKER TYPES ####
# Choose which removal set
sfon_remove <- read.csv(file = "02_data/sfon_nnxnp_markers.csv", header = F) # to remove nnxnp markers
sfon_remove <- as.numeric(sfon_remove)
str(sfon_remove)

length(sfon_remove) # markers to remove
sfon_limited <- drop.markers(sfon, sfon_remove) # drop markers # note that some may already be removed due to seg distortion filter
rm(sfon_remove)

summary(sfon_limited) # check markers
plot.map(sfon_limited, alternate.chrid = T,
         main = "", xlab = "Linkage Group")
text(x = 35, y = 260, labels = paste("n =", sum(nmar(sfon_limited)), "markers"))

# note that re-estimating recombination fraction after dropping those markers doesn't change map positions or total lengths etc.

# sfon <- est.rf(sfon) #doesn't seem to be necessary for crossover determination or for cM position

##export subset data
save.image(file = "02_data/sfon_01_output_subset_only_efxeg_and_hkxhk.RData") # to use for heterochiasmy