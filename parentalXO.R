# An adaptation of the plotGeno() function in order to obtain the parental crossover locations

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
      recalc.chr.length <<- c(recalc.chr.length, max(map))
      
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

# produces an object 'mxoloc.per.chr' and 'dxoloc.per.chr' 
# a dataframe with ind(ividual) and loc(ation))
# mxoloc.per.chr
# dxoloc.per.chr


# Testing using a single round (chr 7, ind 1:10):
parentalXO(sfqtl, chr = 7, ind = c(1:10))



#### 1. OBTAIN parentalXO in all chromosomes ####

# User variables
# For running 'all'
indiv <- 1:nind(sfqtl)
chr <- 1:nchr(sfqtl)
chr <- metacentrics.sfqtl <- c(9,18,1,15,4,7,11,5) # performing for only metacentrics
# indiv <- 1:nind(sfqtl) # still in testing phase

# NULL variables
cum.mxoloc.list <- list(NULL)
cum.dxoloc.list <- list(NULL)
recalc.chr.length <- NULL
cum.recalc.chr.length <- list(NULL)

for(i in chr) {
  parentalXO(sfqtl, chr = i, ind = indiv)
    cum.mxoloc.list[[i]] <- mxoloc.per.chr
    cum.dxoloc.list[[i]] <- dxoloc.per.chr
    cum.recalc.chr.length[[i]] <- recalc.chr.length
}

# Outputs
cum.mxoloc.list
cum.dxoloc.list
cum.recalc.chr.length # Note: not sure why I get multiple numbers, but the final one is the one that is plotted

# Because my P1 = Male instead of standard (= female), mxoloc = P1..


#### 2. COUNT crossovers in multiple chromosomes ####
# Uses local extension of distance for detecting double crossovers

# Choose data variable using either cum.mxoloc.list (FATHER) or cum.dxoloc.list (MOTHER)
#data <- cum.dxoloc.list
#data <- cum.mxoloc.list

# none of this is working right, try to subset earlier...
# metacentrics.sfqtl <- c(9,18,1,15,4,7,11,5)
# data <- cum.dxoloc.list[metacentrics.sfqtl] #metacentric, maternal
# data <- cum.mxoloc.list[metacentrics.sfqtl] #metacentric, paternal
# data <- cum.dxoloc.list[-metacentrics.sfqtl] #acrocentric, maternal
# data <- cum.mxoloc.list[-metacentrics.sfqtl] #acrocentric, paternal


# loop variables
lower <- NULL; upper <- NULL ; odd.even <- NULL ; current.piece <- NULL
indiv.nums <- NULL; counter <- 0 ; XO.spot <- NULL
XO.tot.leng <- NULL ; XO.tot.leng <- NULL ; CUMULATIVE.CHR <- NULL

# user variables
chr <- 1:length(data) # for all
#chr <- c(9,18,1,15,4,7,11,5)
dist <- 100 # distance of 100 cM added to each side

# Create a subset piece from the total list per chromosome
for(i in chr) {
  test <- data[[i]] # temporary to keep below
  print(c("segment", test))
  
  # For each chr, loop to count number of XO
  indiv.nums <- unique(test[,1]) # identify the unique sample names
  print(indiv.nums)
  XO.tot.leng <- c(XO.tot.leng,  cum.recalc.chr.length[[(i)]][length(cum.recalc.chr.length[[(i)]])])
  print("LENGTH OF CHR OF INTEREST")
  print(XO.tot.leng) # use this below to find the relative position of the XO across the CHR
  print("THIS ROUND THE CHR IS")
  print(XO.tot.leng[i])
  current.chr.leng <- XO.tot.leng[i]
  
  for(i in indiv.nums) {
    # find row containing the indiv of interest (per loop)
    indiv.row <- which(test[,1] == (i))
    print(c("indiv.row",indiv.row))
    
    # each row of the indiv of interest
    for(loc in indiv.row) {
      print(test[(loc),2])
      lower <- test[(loc),2] - dist
      upper <- test[(loc),2] + dist
      print(c(lower, upper))
      
      # find the rows of the indiv of interest within the range
      test[indiv.row, 2] > lower & test[indiv.row, 2] < upper 
      print( test[indiv.row, 2] > lower & test[indiv.row, 2] < upper )
      
      # How many other XO for this indiv, chr are within range of this loop's XO 
      print(length(test[indiv.row, 2] > lower & test[indiv.row, 2] < upper)%%2)
      
      # Add the odd (1) or even (0) surrounding XO (within range) value to this loop's XO
      #odd.even <- length(test[indiv.row, 2] > lower & test[indiv.row, 2] < upper)%%2
      odd.even <- table(test[indiv.row, 2] > lower & test[indiv.row, 2] < upper)["TRUE"]%%2
      print(c("odd.even", odd.even))
      
      # With an odd number (1), add a XO to the counter; if even (0) do not add
      if(odd.even == 0) {
        print("even")
        counter <- counter
      } else {
        print("odd")
        print("THIS CHR TOTAL LENGTH IS")
        print(i)
        print(XO.tot.leng[i])
        print("HEY LOOK HERE, IS THIS THE RIGHT LENGTH?")
        print(current.chr.leng)
        XO.spot <<- c(XO.spot, print(test[(loc),2]))
        CUMULATIVE.CHR <<- c(CUMULATIVE.CHR, current.chr.leng)
        counter <- counter + 1 }
      print(c("counter", counter))
    }
  }
}

counter

par(mfrow=c(2,1))
# maternal
hist(XO.spot/CUMULATIVE.CHR*100, xlab = "Relative position of XO (%)", main = "", ylim = c(0,1000))
text(x = 15, y = 900, labels = paste("maternal: n =", length(XO.spot), "XO"))

# paternal
hist(XO.spot/CUMULATIVE.CHR*100, xlab = "Relative position of XO (%)", main = "", ylim = c(0,1000))
text(x = 15, y = 900, labels = paste("paternal: n =", length(XO.spot), "XO"))
