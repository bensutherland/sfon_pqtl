# Count crossovers controlling for likely false double crossovers
# B. Sutherland, labo Bernatchez 2016-10-17
# v0.1
# NOTE: because my P1 = Male instead of standard (= female), mxoloc = P1..

# rm(list=ls())
library(qtl)


##### 0. Create formula #####
# Set NULL for parentalXO
recalc.chr.length <- NULL

# An adaptation of the plotGeno() (R/qtl) function in order to obtain the parental crossover locations
# This produces an object 'mxoloc.per.chr' and 'dxoloc.per.chr', a df w/ ind and loc

parentalXO <- function (x, chr, ind, include.xo = TRUE, horizontal = TRUE, 
                        cutoff = 4, min.sep = 2, cex = 1.2, ...) 
{
  cross <- x
  if (!any(class(cross) == "cross")) 
    stop("Input should have class \"cross\".")
  if (missing(chr)) 
    chr <- names(cross$geno)[1]
  cross <- subset(cross, chr = chr)
  if (nchr(cross) > 1) 
    cross <- subset(cross, chr = names(cross$geno)[1])
  if (!missing(ind)) {
    if (is.null(getid(cross))) 
      cross$pheno$id <- 1:nind(cross)
    if (!is.logical(ind)) 
      ind <- unique(ind)
    cross <- subset(cross, ind = ind)
  }
  id <- getid(cross)
  if (is.null(id)) 
    id <- 1:nind(cross)
  use.id <- TRUE
  type <- class(cross)[1]
  old.las <- par("las")
  on.exit(par(las = old.las))
  par(las = 1)
  if (!("errorlod" %in% names(cross$geno[[1]]))) {
    warning("First running calc.errorlod.")
    cross <- calc.errorlod(cross, error.prob = 0.01)
  }
  errors <- matrix(0, ncol = ncol(cross$geno[[1]]$data), nrow = nrow(cross$geno[[1]]$data))
  dimnames(errors) <- dimnames(cross$geno[[1]]$data)
  top <- top.errorlod(cross, names(cross$geno)[1], cutoff, 
                      FALSE)
  if (length(top) > 0) 
    for (i in 1:nrow(top)) errors[match(top[i, 2], id), as.character(top[i, 
                                                                         3])] <- 1
  map <- cross$geno[[1]]$map
  if (is.matrix(map)) 
    map <- map[1, ]
  L <- diff(range(map))
  min.d <- L * min.sep/100
  d <- diff(map)
  d[d < min.d] <- min.d
  map <- cumsum(c(0, d))
  cross$geno[[1]]$map <- map
  n.ind <- nrow(errors)
  color <- c("white", "gray60", "black", "green", "orange", 
             "red")
  data <- cross$geno[[1]]$data
  chrtype <- class(cross$geno[[1]])
  if (chrtype == "X" && (type == "f2" || type == "bc")) 
    data <- reviseXdata(type, sexpgm = getsex(cross), geno = data, 
                        cross.attr = attributes(cross), force = TRUE)
  if (include.xo) {
    if (type != "4way") {
      xoloc <- locateXO(cross)
      xoloc <- data.frame(ind = rep(1:length(xoloc), sapply(xoloc, 
                                                            length)), loc = unlist(xoloc), stringsAsFactors = TRUE)
    }
    else {
      mcross <- dcross <- cross
      class(mcross)[1] <- class(dcross)[1] <- "bc"
      mcross$geno[[1]]$data[!is.na(data) & data == 1 | 
                              data == 3 | data == 5] <- 1
      mcross$geno[[1]]$data[!is.na(data) & data == 2 | 
                              data == 4 | data == 6] <- 2
      mcross$geno[[1]]$data[!is.na(data) & data == 7 | 
                              data == 8 | data == 9 | data == 10] <- NA
      dcross$geno[[1]]$data[!is.na(data) & data == 1 | 
                              data == 2 | data == 7] <- 1
      dcross$geno[[1]]$data[!is.na(data) & data == 3 | 
                              data == 4 | data == 8] <- 2
      dcross$geno[[1]]$data[!is.na(data) & data == 5 | 
                              data == 6 | data == 9 | data == 10] <- NA
      mxoloc <- locateXO(mcross)
      mxoloc <- data.frame(ind = rep(1:length(mxoloc), 
                                     sapply(mxoloc, length)), loc = unlist(mxoloc), 
                           stringsAsFactors = TRUE)
      dxoloc <- locateXO(dcross)
      dxoloc <- data.frame(ind = rep(1:length(dxoloc), 
                                     sapply(dxoloc, length)), loc = unlist(dxoloc), 
                           stringsAsFactors = TRUE)
      
      ## BEN EDIT START ##
      print("mxoloc")
      print(mxoloc)
      mxoloc.per.chr <<- mxoloc
      
      print("dxoloc")
      print(dxoloc)
      dxoloc.per.chr <<- dxoloc
      
      print("NEXT IS MAP LENGTH")
      print(max(map))
      #recalc.chr.length <<- c(recalc.chr.length, max(map)) # ORIGINAL (WORKS)
      recalc.chr.length <<- max(map) # NEW ATTEMPT
      
      ## BEN EDIT FINISH ##
    }
  }
  args <- list(...)
  if ("main" %in% names(args)) 
    themain <- args$main
  else themain <- paste("Chromosome", names(cross$geno)[1])
  if ("xlim" %in% names(args)) 
    thexlim <- args$xlim
  else thexlim <- NULL
  if ("ylim" %in% names(args)) 
    theylim <- args$ylim
  else theylim <- NULL
  if (type == "4way") {
    jit <- 0.15
    mdata <- data
    ddata <- data
    mdata[!is.na(data) & (data == 1 | data == 3 | data == 
                            5)] <- 1
    mdata[!is.na(data) & (data == 2 | data == 4 | data == 
                            6)] <- 2
    mdata[!is.na(data) & (data == 7 | data == 8)] <- NA
    ddata[!is.na(data) & (data == 1 | data == 2 | data == 
                            7)] <- 1
    ddata[!is.na(data) & (data == 3 | data == 4 | data == 
                            8)] <- 2
    ddata[!is.na(data) & (data == 5 | data == 6)] <- NA
    if (horizontal) {
      if (is.null(thexlim)) 
        thexlim <- c(0, max(map))
      if (is.null(theylim)) 
        theylim <- c(n.ind + 1, 0)
      plot(0, 0, type = "n", xlab = "Location (cM)", ylab = "Individual", 
           main = themain, ylim = theylim, xlim = thexlim, 
           yaxt = "n", yaxs = "i")
      segments(0, 1:n.ind - jit, max(map), 1:n.ind - jit)
      segments(0, 1:n.ind + jit, max(map), 1:n.ind + jit)
      if (use.id) 
        axis(side = 2, at = 1:n.ind, labels = id)
      else axis(side = 2, at = 1:n.ind)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 1] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind - jit, pch = 21, col = "black", bg = color[1], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 2] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind - jit, pch = 21, col = "black", bg = color[3], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 9] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind - jit, pch = 21, col = "black", bg = color[4], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 10] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind - jit, pch = 21, col = "black", bg = color[5], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 1] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind + jit, pch = 21, col = "black", bg = color[1], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 2] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind + jit, pch = 21, col = "black", bg = color[3], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 9] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind + jit, pch = 21, col = "black", bg = color[4], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 10] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind + jit, pch = 21, col = "black", bg = color[5], 
             cex = cex)
      u <- par("usr")
      segments(map, u[3], map, u[3] - 1/2)
      segments(map, u[4], map, u[4] + 1/2)
      if (any(errors != 0)) {
        ind <- rep(1:n.ind, length(map))
        ind[errors != 1] <- NA
        points(x, ind - jit, pch = 0, col = color[6], 
               cex = cex + 0.4, lwd = 2)
        points(x, ind + jit, pch = 0, col = color[6], 
               cex = cex + 0.4, lwd = 2)
      }
      if (include.xo) {
        points(mxoloc$loc, mxoloc$ind - jit, pch = 4, 
               col = "blue", lwd = 2)
        points(dxoloc$loc, dxoloc$ind + jit, pch = 4, 
               col = "blue", lwd = 2)
      }
    }
    else {
      if (is.null(theylim)) 
        theylim <- c(max(map), 0)
      if (is.null(thexlim)) 
        thexlim <- c(0, n.ind + 1)
      plot(0, 0, type = "n", ylab = "Location (cM)", xlab = "Individual", 
           main = themain, xlim = thexlim, ylim = theylim, 
           xaxt = "n", xaxs = "i")
      segments(1:n.ind - jit, 0, 1:n.ind - jit, max(map))
      segments(1:n.ind + jit, 0, 1:n.ind + jit, max(map))
      if (use.id) 
        axis(side = 1, at = 1:n.ind, labels = id)
      else axis(side = 1, at = 1:n.ind)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 1] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind - jit, y, pch = 21, col = "black", bg = color[1], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 2] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind - jit, y, pch = 21, col = "black", bg = color[3], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 9] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind - jit, y, pch = 21, col = "black", bg = color[4], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(mdata)] <- NA
      ind <- tind
      ind[!is.na(mdata) & mdata != 10] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind - jit, y, pch = 21, col = "black", bg = color[5], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 1] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind + jit, y, pch = 21, col = "black", bg = color[1], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 2] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind + jit, y, pch = 21, col = "black", bg = color[3], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 9] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind + jit, y, pch = 21, col = "black", bg = color[4], 
             cex = cex)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(ddata)] <- NA
      ind <- tind
      ind[!is.na(ddata) & ddata != 10] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind + jit, y, pch = 21, col = "black", bg = color[5], 
             cex = cex)
      u <- par("usr")
      segments(u[1], map, (u[1] + 1)/2, map)
      segments(u[2], map, (n.ind + u[2])/2, map)
      if (any(errors != 0)) {
        ind <- rep(1:n.ind, length(map))
        ind[errors != 1] <- NA
        points(ind - jit, y, pch = 0, col = color[6], 
               cex = cex + 0.4, lwd = 2)
        points(ind + jit, y, pch = 0, col = color[6], 
               cex = cex + 0.4, lwd = 2)
      }
      if (include.xo) {
        points(mxoloc$ind - jit, mxoloc$loc, pch = 4, 
               col = "blue", lwd = 2)
        points(dxoloc$ind + jit, dxoloc$loc, pch = 4, 
               col = "blue", lwd = 2)
      }
    }
  }
  else {
    if (horizontal) {
      if (is.null(thexlim)) 
        thexlim <- c(0, max(map))
      if (is.null(theylim)) 
        theylim <- c(n.ind + 0.5, 0.5)
      plot(0, 0, type = "n", xlab = "Location (cM)", ylab = "Individual", 
           main = themain, ylim = theylim, xlim = thexlim, 
           yaxt = "n")
      segments(0, 1:n.ind, max(map), 1:n.ind)
      if (use.id) 
        axis(side = 2, at = 1:n.ind, labels = id)
      else axis(side = 2)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(data)] <- NA
      ind <- tind
      ind[!is.na(data) & data != 1] <- NA
      x <- rep(map, rep(n.ind, length(map)))
      points(x, ind, pch = 21, col = "black", bg = color[1], 
             cex = cex)
      ind <- tind
      ind[!is.na(data) & data != 2] <- NA
      if (type == "f2" || (type == "bc" && chrtype == "X")) 
        points(x, ind, pch = 21, col = "black", bg = color[2], 
               cex = cex)
      else points(x, ind, pch = 21, col = "black", bg = color[3], 
                  cex = cex)
      if (type == "f2" || (type == "bc" && chrtype == "X")) {
        ind <- tind
        ind[!is.na(data) & data != 3] <- NA
        points(x, ind, pch = 21, col = "black", bg = color[3], 
               cex = cex)
      }
      if (type == "f2") {
        ind <- tind
        ind[!is.na(data) & data != 4] <- NA
        points(x, ind, pch = 21, col = "black", bg = color[4], 
               cex = cex)
        ind <- tind
        ind[!is.na(data) & data != 5] <- NA
        points(x, ind, pch = 21, col = "black", bg = color[5], 
               cex = cex)
      }
      u <- par("usr")
      segments(map, u[3], map, u[3] - 1/2)
      segments(map, u[4], map, u[4] + 1/2)
      if (any(errors != 0)) {
        ind <- rep(1:n.ind, length(map))
        ind[errors != 1] <- NA
        points(x, ind, pch = 0, col = color[6], cex = cex + 
                 0.4, lwd = 2)
      }
      if (include.xo) 
        points(xoloc$loc, xoloc$ind, pch = 4, col = "blue", 
               lwd = 2)
    }
    else {
      if (is.null(theylim)) 
        theylim <- c(max(map), 0)
      if (is.null(thexlim)) 
        thexlim <- c(0.5, n.ind + 0.5)
      plot(0, 0, type = "n", ylab = "Location (cM)", xlab = "Individual", 
           main = themain, xlim = thexlim, ylim = theylim, 
           xaxt = "n")
      segments(1:n.ind, 0, 1:n.ind, max(map))
      if (use.id) 
        axis(side = 1, at = 1:n.ind, labels = id)
      else axis(side = 1)
      tind <- rep(1:n.ind, length(map))
      tind[is.na(data)] <- NA
      ind <- tind
      ind[!is.na(data) & data != 1] <- NA
      y <- rep(map, rep(n.ind, length(map)))
      points(ind, y, pch = 21, col = "black", bg = "white", 
             cex = cex)
      ind <- tind
      ind[!is.na(data) & data != 2] <- NA
      if (type == "f2" || (type == "bc" && chrtype == "X")) 
        points(ind, y, pch = 21, col = "black", bg = color[2], 
               cex = cex)
      else points(ind, y, pch = 21, col = "black", bg = color[3], 
                  cex = cex)
      if (type == "f2" || (type == "bc" && chrtype == "X")) {
        ind <- tind
        ind[!is.na(data) & data != 3] <- NA
        points(ind, y, pch = 21, col = "black", bg = color[3], 
               cex = cex)
      }
      if (type == "f2") {
        ind <- tind
        ind[!is.na(data) & data != 4] <- NA
        points(ind, y, pch = 21, col = "black", bg = color[4], 
               cex = cex)
        ind <- tind
        ind[!is.na(data) & data != 5] <- NA
        points(ind, y, pch = 21, col = "black", bg = color[5], 
               cex = cex)
      }
      u <- par("usr")
      segments(u[1], map, (u[1] + 1)/2, map)
      segments(u[2], map, (n.ind + u[2])/2, map)
      if (any(errors != 0)) {
        ind <- rep(1:n.ind, length(map))
        ind[errors != 1] <- NA
        points(ind, y, pch = 0, col = color[6], cex = cex + 
                 0.4, lwd = 2)
      }
      if (include.xo) 
        points(xoloc$ind, xoloc$loc, pch = 4, col = "blue", 
               lwd = 2)
    }
  }
  invisible()
}


###### Load Data ######
# set working directory
setwd("~/Documents/bernatchez/01_Sfon_projects/03_Sfon_pQTL/sfon_pqtl")

# Load Part 1 results:
# load("02_data/sfon_01_output.RData")
load("02_data/sfon_01_output_subset_only_efxeg_and_hkxhk.RData")
sfon_limited

# Identify male and female individuals (plotGeno uses 'ind' arg)
ind.males = c(sfon$pheno$sex=="M")
ind.females = c(sfon$pheno$sex=="F")

# Identify autosomes, sex chr, metacentrics and acrocentrics (plotGeno uses 'chr' arg)
metacentrics <- c(1:8)
acrocentrics <- c(9:34,36:42)
sex.chrom <- 35
autosomes <- c(1:34,36:42)

#### OBTAIN parentalXO in sets ####

# create an object with directions of which chromosomes to obtain
cross <- NULL; indiv <- NULL; chr = NULL
outing <- NULL; outing <- list(); name.of.sets <- NULL
collect.me <- NULL; collect.me <- list()
index <- NULL
sets <- NULL; sets <- list()

# Choose either a standard run to average over metacentrics and acrocentrics
# or, choose the special version that runs all chromosomes individually (see below)

# Standard
sets[[1]] <- metacentrics
sets[[2]] <- acrocentrics
sets[[3]] <- sex.chrom
name.of.sets <- c("metacentrics","acrocentrics","sex.chrom") # temporary needs to be fix

# Special (all chromosomes individually)
# sets[[1]] <- 1; sets[[2]] <- 2 ; sets[[3]] <- 3; sets[[4]] <- 4 ; sets[[5]] <- 5; sets[[6]] <- 6; sets[[7]] <- 7
# sets[[8]] <- 8; sets[[9]] <- 9; sets[[10]] <- 10; sets[[11]] <- 11; sets[[12]] <- 12; sets[[13]] <- 13
# sets[[14]] <- 14; sets[[15]] <- 15; sets[[16]] <- 16; sets[[17]] <- 17; sets[[18]] <- 18
# sets[[19]] <- 19; sets[[20]] <- 20; sets[[21]] <- 21; sets[[22]] <- 22; sets[[23]] <- 23
# sets[[24]] <- 24; sets[[25]] <- 25; sets[[26]] <- 26; sets[[27]] <- 27; sets[[28]] <- 28
# sets[[29]] <- 29; sets[[30]] <- 30; sets[[31]] <- 31; sets[[32]] <- 32; sets[[33]] <- 33
# sets[[34]] <- 34; sets[[35]] <- 35; sets[[36]] <- 36; sets[[37]] <- 37; sets[[38]] <- 38
# sets[[39]] <- 39; sets[[40]] <- 40; sets[[41]] <- 41; sets[[42]] <- 42
# name.of.sets <- seq(1:42)

# Loop to subset 
for(i in 1:length(sets)) {
  set.name <- NULL
  set.name <- name.of.sets[i]
  print(set.name) #provide a name for the set
  
  #subset data
  cross <- subset(sfon_limited, chr = sets[[i]])
  print(nchr(cross))
  indiv <- 1:nind(cross) # for full set
  #indiv <- 1:10
  chr <- as.numeric(names(cross$geno))
  
  # Obtain parentalXO in selected set
  # set NULL
  cum.mxoloc.list <- list(NULL) ; cum.dxoloc.list <- list(NULL) ; recalc.chr.length <- NULL 
  cum.recalc.chr.length <- NULL; name <- NULL
  
  for(c in chr) {
    parentalXO(cross, chr = c, ind = indiv)
    
    # new method that instead of using i uses name
    name <- paste("chr",c,sep="") # BETTER
    print(name)
    
    cum.mxoloc.list[[name]] <- mxoloc.per.chr
    cum.dxoloc.list[[name]] <- dxoloc.per.chr
    cum.recalc.chr.length[[name]] <- recalc.chr.length
  }
  
  str(cum.mxoloc.list)
  str(cum.dxoloc.list)
  str(cum.recalc.chr.length)
  
  #### Count and correct for Double XO #####
  # user variables
  distance <- 50 # distance to be screened on each side for crossovers
  #distance <- 0.1 # if want to evaluate what happens if no double crossovers are removed
  
  # Choose data variable using either cum.dxoloc.list (here: MOTHER) or cum.mxoloc.list (here: FATHER)
  
  both.data <- NULL
  both.data <- list()
  both.data[[1]] <- cum.dxoloc.list
  both.data[[2]] <- cum.mxoloc.list
  names(both.data) <- c("cum.dxoloc.list","cum.mxoloc.list")
  
  for(d in 1:length(both.data)){
    data <- both.data[[d]]
    print(data)
    name.of.data <- NULL
    name.of.data <- names(both.data)[d]
    print(name.of.data)

      # Set up to obtain the crossovers
      # NULL variables
      lower <- NULL; upper <- NULL ; odd.even <- NULL ; current.piece <- NULL
      indiv.nums <- NULL; counter <- 0 ; XO.spot <- NULL
      XO.tot.leng <- NULL ; CUMULATIVE.CHR <- NULL
      per.chromosome.XO <- NULL
      name <- NULL
      even.counter <- 0
      
      # Create a subset piece from the total list per chromosome
      for(n in chr) {
        
        # new method that instead of using i uses name
        name <- paste("chr",n,sep="")
        print(name)
        
        # EXPERIMENTAL
        test <- data[[name]] # take out data from one chromosome
        print(c("*Treating chromosome", n), quote = F)
        print(test)
        per.chromosome.counter <- 0
        
        # For each chromosome, count the number of XO per individual
        indiv.nums <- unique(test[,1]) # identify the unique sample names in 'test'
        print(c("**Treating Sample:", indiv.nums), quote=F)
        
        # For each chromosome, record the total length
        XO.tot.leng <- c(XO.tot.leng,  cum.recalc.chr.length[name])
        print(c("***XO.tot.len", XO.tot.leng), quote=F)  
        
        # Extract the length of this chromosome
        current.chr.leng <- NULL
        current.chr.leng <- XO.tot.leng[name]
        print(c("THIS ROUND THE CHR IS", XO.tot.leng[name]))
        
        # Per chromosome, check each unique indiv
        for(j in indiv.nums) {
          indiv.row <- which(test[,1] == (j)) # find row(s) for the indiv of interest (per loop)
          print(c("THIS IS ****", j))
          print(c("indiv.row",indiv.row))
          
          # Per unique indiv, check each XO (i.e. how many within range?)
          for(loc in indiv.row) {
            
            print("DEFINING RANGE")
            print(test[(loc),2])
            lower <- test[(loc),2] - distance
            upper <- test[(loc),2] + distance
            print(c(lower, upper))
            
            # TRUE for each XO in range
            print("FOR EACH INDIVIDUAL IN THIS CHR, score TRUE for each XO within range")
            test[indiv.row, 2] > lower & test[indiv.row, 2] < upper 
            print( test[indiv.row, 2] > lower & test[indiv.row, 2] < upper )
            
            # Number of XO within the range
            print(c(length(test[indiv.row, 2] > lower & test[indiv.row, 2] < upper), "XOs on chr")) 
            
            # 0 if EVEN, 1 if ODD
            print("**** 0 if EVEN, 1 if ODD ****")
            print(length(test[indiv.row, 2] > lower & test[indiv.row, 2] < upper)%%2) 
            
            # True to add 1 to counter
            odd.even <- table(test[indiv.row, 2] > lower & test[indiv.row, 2] < upper)["TRUE"]%%2
            print("IS THIS XO ODD OR EVEN")
            print(odd.even)
            
            # If odd number (1), add a XO to the counter; if even (0) do not add
            if(odd.even == 0) {
              print("even")
              counter <- counter
              even.counter <- even.counter + 1
            } else {
              print("odd")
              
              # Collect the location of this crossover
              XO.spot <<- c(XO.spot, print(test[(loc),2]))
              
              # Collect the total chromosome length for this crossover
              CUMULATIVE.CHR <<- c(CUMULATIVE.CHR, current.chr.leng)
              #print(c("CUMULATIVE.CHR", CUMULATIVE.CHR))
              
              # Add one to counter
              counter <- counter + 1 
              
              # Add one to per.chromosome.counter
              per.chromosome.counter <- per.chromosome.counter + 1
            }
            print(c("counter", counter))
          }
        }
        print(c("per.chromosome.counter"))
        print(c(per.chromosome.counter, "chr", i ))
        per.chromosome.XO[name] <- per.chromosome.counter
      }
      
      counter
      even.counter
  
  # USE THE FOLLOWING VARIABLES TO CREATE OBJECTS: 
  name.of.data # dxoloc or mxoloc? (stillneed to make for loop)
  set.name # sets piece
  
  # obtain this information
  collect.me.temp <- XO.spot/CUMULATIVE.CHR
  index <- paste(set.name, name.of.data, sep="_")
  collect.me[[index]] <- collect.me.temp
  
}
}



# This should have all you need
str(collect.me)

names(collect.me)

# Plot either standard or Special that runs each individual chr separately

#Standard
par(mfrow=c(2,3), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
order <- c(1,3,5,2,4,6)
1:length(names(collect.me))

# #Special (individual chromosomes)
# par(mfrow=c(6,7), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
# order <- c(seq(from = 1, to = 84, by = 2)
#            , seq(from = 2, to = 84, by = 2))
# length(order)

# Plot
maximum <- 0.10
#maximum <- 500 # when freq = T

for(p in order){
  chr.set <- collect.me[[p]]
  hist(chr.set*100, xlab = "Relative position of XO (%)", main = ""
       , xlim = c(0,100), las = 1
       , ylim = c(0,maximum)
       # or
       , freq = FALSE# to display as percent
       , breaks = 10
      )
  text(x = 50, y = maximum-0.01, labels = paste("avg/chr ="
                                                , round(
                                                  length(chr.set)/length(unique(names(chr.set)))
                                                  ,1)))
  text(x= 50, y = maximum, labels = names(collect.me)[p])
}

# save as 7.3 * 3.8 in portrait
#or
# save as 7 * 4.5 for six panel


# save as 14 by 14 when all indiv





#### HERE BE DRAGONS #####
#### Summary Statistics ####
#set variables
#sfqtl
# metacentrics <- c(9,18,1,15,4,7,11,5) #only metacentrics
# acrocentrics <- c(2:3,6,8,10,12:14,16:17,19:42) #acrocentrics
# sex.chrom <- 22
# autosomes <- c(1:21,23:42)
# 
# 
# # sum, average, and 
# print(c("sum"))
# sum(per.chromosome.XO[c(metacentrics)])
# print(c("average"))
# mean(per.chromosome.XO[c(metacentrics)])
# print(c("sd"))
# sd(per.chromosome.XO[c(metacentrics)])
# 
# print(c("sum"))
# sum(per.chromosome.XO[c(acrocentrics)])
# print(c("average"))
# mean(per.chromosome.XO[c(acrocentrics)])
# print(c("sd"))
# sd(per.chromosome.XO[c(acrocentrics)])
# 
# print(c("sum"))
# sum(per.chromosome.XO[c(sex.chrom)])
# print(c("average"))
# mean(per.chromosome.XO[c(sex.chrom)])
# print(c("sd"))
# sd(per.chromosome.XO[c(sex.chrom)])
# 
# print(c("sum"))
# sum(per.chromosome.XO[c(autosomes)])
# print(c("average"))
# mean(per.chromosome.XO[c(autosomes)])
# print(c("sd"))
# sd(per.chromosome.XO[c(autosomes)])
# 
# 
# ## attempt to automate summary statistics
# # # summary statistics
# # sum(per.chromosome.XO)[metacentrics]
# # 
# # # set up a vector for your conformation
# # chr.conformation <- NULL
# # chr.conformation[metacentrics] <- "metacentric"
# # chr.conformation[acrocentrics] <- "acrocentric"
# # 
# # # make a dataframe
# # per.chromosome.XO.df <- as.data.frame(cbind(per.chromosome.XO, chr.conformation), row.names = 1:42)
# # 
# # # perform some tests
# # mean(per.chromosome.XO.df[,2])