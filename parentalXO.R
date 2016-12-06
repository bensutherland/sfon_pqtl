# Count crossovers controlling for likely false double crossovers
# B. Sutherland, labo Bernatchez 2016-10-17
# v0.1

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
acrocentrics <- c(9:42)
sex.chrom <- 35
autosomes <- c(1:34,36:42)

# Create cross types with the desired chromosomes
#sfon_limited_metacentrics <- subset(sfon_limited, chr = metacentrics)
#sfon_limited_acrocentrics <- subset(sfon_limited, chr = acrocentrics)
#sfon_limited_acrocentrics <- subset(sfon_limited, chr = sex.chrom)
#sfon_limited_acrocentrics <- subset(sfon_limited, chr = autosomes)


##### Set Data #####
cross <- sfon_limited # for efxeg and hkxhk markers only
#cross <- sfon_limited_metacentrics # for metacentrics
#cross <- sfon_limited_acrocentrics # for acrocentrics
#cross <- sfon_limited_sex_chrom
#cross <- sfon_limited_autosomes


# Select which individuals to include
indiv <- 1:nind(cross) # all individuals
#chr <- as.numeric(names(cross$geno))
chr <- 10:20
#chr <- 41

# first need to calculate errorlod
# sfon_limited <- calc.errorlod(cross = cross)
# note that this will change the lengths in the map

# Test to make sure formula works using a single round (chr 7, ind 1:10):
parentalXO(cross, chr = 7, ind = c(1:10))


##### TROUBLESHOOTING STARTS HERE ###
# # mini tests
# chr <- 2
# indiv <- 1:20


#### 1. OBTAIN parentalXO in all chromosomes ####
# NULL variables
cum.mxoloc.list <- list(NULL) ; cum.dxoloc.list <- list(NULL) ; recalc.chr.length <- NULL ; cum.recalc.chr.length <- NULL
name <- NULL

for(i in chr) {
    parentalXO(cross, chr = i, ind = indiv)
    
    # new method that instead of using i uses name
    name <- paste("chr",i,sep="")
    print(name)
    
    cum.mxoloc.list[[name]] <- mxoloc.per.chr
    cum.dxoloc.list[[name]] <- dxoloc.per.chr
    cum.recalc.chr.length[[name]] <- recalc.chr.length
}

# Outputs
str(cum.mxoloc.list)
str(cum.dxoloc.list)

cum.recalc.chr.length

# NOTE: because my P1 = Male instead of standard (= female), mxoloc = P1..



#### 2. COUNT crossovers in multiple chromosomes ####
# Uses local extension of distance for detecting double crossovers

# user variables
distance <- 50 # distance to be screened on each side for crossovers
#chr <- 1:length(data) # for all

# Choose data variable using either cum.mxoloc.list (here: FATHER) or cum.dxoloc.list (here: MOTHER)
data <- cum.dxoloc.list #P2
#data <- cum.mxoloc.list #P1


# NULL variables
lower <- NULL; upper <- NULL ; odd.even <- NULL ; current.piece <- NULL
indiv.nums <- NULL; counter <- 0 ; XO.spot <- NULL
XO.tot.leng <- NULL ; CUMULATIVE.CHR <- NULL
per.chromosome.XO <- NULL
name <- NULL

# Create a subset piece from the total list per chromosome
for(i in chr) {
  
  # new method that instead of using i uses name
  name <- paste("chr",i,sep="")
  print(name)
  
  # EXPERIMENTAL
  test <- data[[name]] # take out data from one chromosome
  
  #test <- data[[i]] # take out data from one chromosome
  print(c("*Treating chromosome", i), quote = F)
  print(test)
  per.chromosome.counter <- 0
  
  # For each chromosome, count the number of XO per individual
  indiv.nums <- unique(test[,1]) # identify the unique sample names in 'test'
  print(c("**Treating Sample:", indiv.nums), quote=F)
  
  # For each chromosome, record the total length
  XO.tot.leng <- c(XO.tot.leng,  cum.recalc.chr.length[name])
  print(c("***XO.tot.len", XO.tot.leng), quote=F)  
  
  # extract this chromosome loop's chromosome length
  current.chr.leng <- NULL
  current.chr.leng <- XO.tot.leng[name]
  print(c("THIS ROUND THE CHR IS", XO.tot.leng[name]))
  
  # NULL this chromosome's counter
  #per.chromosome.counter <- NULL
  ### IS THIS OK?
    
    # Per chromosome, check each unique indiv
    for(j in indiv.nums) {
      indiv.row <- which(test[,1] == (j)) # find row(s) for the indiv of interest (per loop)
      print(c("THIS IS ****", j))
      print(c("indiv.row",indiv.row))
      
      # Per unique indiv, check each XO
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


#### 3. PLOT the male and female graphs #####
par(mfrow=c(2,1), mar= c(3,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
# maternal
hist(XO.spot/CUMULATIVE.CHR*100, xlab = "Relative position of XO (%)", main = "", ylim = c(0,1000), las = 1 )
text(x = 15, y = 900, labels = paste("maternal: n =", length(XO.spot), "XO"))

# paternal
hist(XO.spot/CUMULATIVE.CHR*100, xlab = "Relative position of XO (%)", main = "", ylim = c(0,1000), las = 1)
text(x = 15, y = 900, labels = paste("paternal: n =", length(XO.spot), "XO"))

# save as 7.3 * 3.8 in portrait




#### Summary Statistics ####
#set variables
#sfqtl
metacentrics <- c(9,18,1,15,4,7,11,5) #only metacentrics
acrocentrics <- c(2:3,6,8,10,12:14,16:17,19:42) #acrocentrics
sex.chrom <- 22
autosomes <- c(1:21,23:42)


# sum, average, and 
print(c("sum"))
sum(per.chromosome.XO[c(metacentrics)])
print(c("average"))
mean(per.chromosome.XO[c(metacentrics)])
print(c("sd"))
sd(per.chromosome.XO[c(metacentrics)])

print(c("sum"))
sum(per.chromosome.XO[c(acrocentrics)])
print(c("average"))
mean(per.chromosome.XO[c(acrocentrics)])
print(c("sd"))
sd(per.chromosome.XO[c(acrocentrics)])

print(c("sum"))
sum(per.chromosome.XO[c(sex.chrom)])
print(c("average"))
mean(per.chromosome.XO[c(sex.chrom)])
print(c("sd"))
sd(per.chromosome.XO[c(sex.chrom)])

print(c("sum"))
sum(per.chromosome.XO[c(autosomes)])
print(c("average"))
mean(per.chromosome.XO[c(autosomes)])
print(c("sd"))
sd(per.chromosome.XO[c(autosomes)])



## attempt to automate summary statistics
# # summary statistics
# sum(per.chromosome.XO)[metacentrics]
# 
# # set up a vector for your conformation
# chr.conformation <- NULL
# chr.conformation[metacentrics] <- "metacentric"
# chr.conformation[acrocentrics] <- "acrocentric"
# 
# # make a dataframe
# per.chromosome.XO.df <- as.data.frame(cbind(per.chromosome.XO, chr.conformation), row.names = 1:42)
# 
# # perform some tests
# mean(per.chromosome.XO.df[,2])



