nna <- function(x) x[!is.na(x)]

# read data from shapefile (SLOW for high resolution)
read.naturalearth <- function(infile){
  require(maptools)
  PW <- readShapePoly(infile)
  ncountries <- length(PW@polygons)
  ngon <- rep(NA,ncountries)
  ndata <- 0
### how many polygons/country
  for(i in 1:ncountries){
    ngon[i] <- length(PW@polygons[[i]]@Polygons)
  }
  totgon <- sum(ngon)
  lLen <- rep(0,totgon)
  k <- 1
### how many data points / polygon
  for(i in 1:ncountries){
    for (j in 1:ngon[i]){
      lLen[k] <-  dim(PW@polygons[[i]]@Polygons[[j]]@coords)[1]
      k <- k+1
    }
  }
  ndata <- sum(lLen) + totgon  # add # polygons for NA seperators!
  lE <- cumsum(lLen+1)-1
  lB <- lE - lLen + 1
  wdata <- data.frame(x=rep(NA,ndata),y=rep(NA,ndata))
### in v2.0.0: names are in lower case
### but v3.1.0: UPPER CASE
  if(is.element("NAME",names(PW@data))) {
    cnames <- as.vector(PW@data$NAME_LONG)
  } else if(is.element("name",names(PW@data))) {
    cnames <- as.vector(PW@data$name_long)
  } else stop("Can't find name field.")

  gon.names <- character(totgon)
  k <- 1
  for(i in 1:ncountries){
    for (j in 1:ngon[i]){
      if(ngon[i]==1) gon.names[k] <- cnames[i]
      else gon.names[k] <- paste(cnames[i],":",j,sep="")
      cc <- data.frame(PW@polygons[[i]]@Polygons[[j]]@coords)
#      cat(i,j,k,lB[k],lE[k],"\n")
      wdata$x[lB[k]:lE[k]] <- cc[,1]
      wdata$y[lB[k]:lE[k]] <- cc[,2]
      k <- k+1
    }
  }
#  wdata <- round(wdata*pi/180*10^6)
# create list of default polyon names from country names list
  gLen <- rep(1,totgon)
  gE <- 2*(1:totgon)
  gB <- gE - 1
  gdata <- rep(NA,2*totgon)
  cat("There are ",totgon,"polygons.\n")
  gdata[seq(1,2*totgon-1,by=2)] <- 1:totgon
  result <- list(line=list(x=wdata$x,y=wdata$y,B=lB,E=lE,length=lLen),
                 gon=list(names=gon.names,B=gB,E=gE,length=gLen,data=gdata),
                 ngon=ngon,ncountries=ncountries,cnames=cnames)
  result
}

# test: get 1 polygon from the .shp file (or all polygons for 1 country)
read.ne <- function(infile,country=1,gon=1){
  require(maptools)
  PW <- readShapePoly(infile)
  print(as.vector(PW@data$name_long)[country])
  if(gon>0) {
    wdata <- data.frame(PW@polygons[[country]]@Polygons[[gon]]@coords)
    names(wdata) <- c("x","y")
  }
  else {
    ngon <- length( PW@polygons[[country]]@Polygons)
    cat("There are",ngon,"gons.\n")
    ww0 <- data.frame(x=NA,y=NA)
    wdata <- NULL
    for (j in 1:ngon){
      cc <- data.frame(PW@polygons[[country]]@Polygons[[j]]@coords)
      names(cc) <- c("x","y")
      wdata <- if(is.null(wdata)) cc else rbind(wdata,ww0,cc)
    }
  }
  wdata
}

# test: common points in 2 polygons
gon.common <- function(g1,g2){
  border <- (g1$x %in% g2$x) & (g1$y %in% g2$y) 
#  print(diff(border))
  if(sum(diff(border)==1) > 1 | sum(diff(border)== -1)>1) warning("Strange border?") 
  data.frame(x=g1$x[border],y=g1$y[border])
}


#2. check valence of lines: 1 is coast, 2 is border, >3 is vertex 
#   COMPLICATION: start/end of line
# loop over all points and see how often they appear in the polygons
# # = 1 : coastal point
# # = 2 : a political border -> appears in 2 polygons
#         OR begin/end point of a polygon
# # > 2 : a vertex : 2 or more borders come together in this point
# useage: ww$valence <- map.valence(ww$data)

# THIS USES C CODE, BUT IT STILL CAN TAKE HALF AN HOUR !
# valence if the values are integers e.g. round(x*10^8)
# somewhat faster than for reals...
map.valence <- function(xy){
  NX=length(xy$x)
  .C("mapvalence",x=as.integer(xy$x),y=as.integer(xy$y),len=as.integer(NX),
                  valence=integer(NX),NAOK=TRUE)$valence
}

# remove duplicate points at given precision
# if the values are integer, you can set prec=0
# INCONSISTENCY: for non-zero precision, you may remove different points
# depending on the direction of the line
# so the border between 2 polygons may becomes "different"
map.clean <- function(xy){
  NX=length(xy$x)
  result <- .C("mapclean",x=as.integer(xy$x),y=as.integer(xy$y),len=as.integer(NX),
               nx=integer(NX),ny=integer(NX),nlen=integer(1),
               NAOK=TRUE)
  data.frame(x=result$nx[1:result$nlen],y=result$ny[1:result$nlen])
}

# split the polygons into lines
# look for *changes in valence* along a set of points.
# These become split points
# BUGS: - a coast line of length 1 is not identified: two consecutive points 
#         are both vertices with valence 2
#         luckily, this only happens a few times, so we correct manually
# ATTENTION: the data MUST start with NA, or the routine will kill R
map.split <- function(xy){
  if(is.null(xy$val)) xy$val <- map.valence(xy)
#  if(xy$val[1] != 0) stop("The data should have NA in the first row!")
  NX <- length(xy$x)
  MAXGON=20*sum(is.na(xy$x)) # average of 20 polylines/gon (it's MUCH less)
  MAXLIN=NX*2 # we expect a slightly longer data set (repeated end points) 
  
  result <- .C("mapsplit",x=as.integer(xy$x),y=as.integer(xy$y),len=as.integer(NX),
               nx=integer(MAXLIN),ny=integer(MAXLIN),nlen=integer(1),
               valence=as.integer(xy$val),gon=integer(MAXGON),ngon=integer(1),NAOK=TRUE)
  
  list(line=data.frame(x=result$nx[1:result$nlen],y=result$ny[1:result$nlen]),
       gon=result$gon[1:result$ngon])
}

# reformat the line/gon data into a canonical map format
make.map <- function(ww){
  gN <- sum(is.na(ww$gon))
  gE <- which(is.na(ww$gon))-1
  gB <- c(1,gE[-gN]+2)
  gL <- gE-gB+1
  newgon <- list(data=ww$gon,names=NULL,
                 ngons=gN,E=gE,B=gB,length=gL)
  lN <- sum(is.na(ww$line$x))
  lE <- which(is.na(ww$line$x))-1
  lB <- c(1,lE[-lN]+2)
  lLen <- lE-lB+1
  
  lR <- rep(0,lN)
  lL <- rep(0,lN)
  for(i in 1:lN) {
    loc <- which(ww$gon==i)
    lL[i] <- which(gB<=loc & gE>=loc)
  }
  newline <- list(x=ww$line$x,y=ww$line$y,L=lL,R=lR,length=lLen,
                  B=lB,E=lE,nlines=lN)
  list(gon=newgon,line=newline)
}

#############
### duplicate lines 
find.dup <- function(loc=1,ww=world,quiet=TRUE){
  len <- ww$line$B[loc]
  xy <- get.line(loc,ww)
  N <- dim(xy)[1]
  candidates <- which(ww$line$length==ww$line$length[loc])
  candidates <- candidates[candidates!=loc]

  duplist <- numeric(0)
  print(candidates)
  for(i in candidates){
    xy2= get.line(i,ww)
    if(!any(xy2!=xy)) {
      if(i!=loc) duplist <- c(duplist,i)
    }
    else if(!any(xy2[N:1,]!=xy)) {
      duplist <- c(duplist,-i)
    }
  }
  duplist
}
# remove duplicate lines
remove.dup <- function(loc=1,ww=world,quiet=TRUE){
  len <- ww$line$B[loc]
  xy <- get.line(loc,ww)
  N <- dim(xy)[1]
  candidates <- which(ww$line$length==ww$line$length[loc])
  candidates <- candidates[candidates!=loc]
  DUP <- NA
  nc <- length(candidates)
  i <- 1
  while(i<=nc & is.na(DUP)){
    cc <- candidates[i]
    xy2= get.line(cc,world)
    if(!any(xy2!=xy)) {
      if (!quiet) print(paste("removing DUPLICATE of ",loc,":",cc))
      DUP <- cc
    }
    else if(!any(xy2[N:1,]!=xy)) {
      if (!quiet) print(paste("removing REVERSE DUPLICATE of ",loc,":",cc))
      DUP <- -cc
    }
    i <- i+1
  }
  if(!is.na(DUP)){
    ww$gon$data[ww$gon$data==abs(DUP)] <- sign(DUP)*loc
    ww$gon$data[ww$gon$data== -abs(DUP)] <- -sign(DUP)*loc
    ww <- remove.line(abs(DUP),ww)
  }
  ww
}

make.LR <- function(ww) {
  lN <- ww$line$nlines
  for(i in 1:lN){
    lL <- which(ww$gon$data==i)
    lR <- which(ww$gon$data== -i)
    if(length(lR)>1 | length(lL)>1) warning(paste("Line",i,"multiple match"))
    ww$line$L[i] <- if(length(lL)>0) which(ww$gon$B<=lL[1] & ww$gon$E>=lL[1]) else 0
    ww$line$R[i] <- if(length(lR)>0) which(ww$gon$B<=lR[1] & ww$gon$E>=lR[1]) else 0
  }
  ww
}

### read a set of lines (e.g. rivers)
### NOT YET FINISHED!!!
read.naturalearth.line <- function(infile,get.data=TRUE){
  require(maptools)
  PW <- readShapeSpatial(infile)
  ncountries <- length(PW@lines)
  ngon <- rep(NA,ncountries)
  ndata <- 0
### how many polygons/country
  for(i in 1:ncountries){
    ngon[i] <- length(PW@lines[[i]]@Lines)
  }
  totgon <- sum(ngon)
  lLen <- rep(0,totgon)
  k <- 1
### how many data points / polygon
  for(i in 1:ncountries){
    for (j in 1:ngon[i]){
      lLen[k] <-  dim(PW@lines[[i]]@Lines[[j]]@coords)[1]
      k <- k+1
    }
  }
  ndata <- sum(lLen) + totgon  # add # polygons for NA seperators!
  lE <- cumsum(lLen+1)-1
  lB <- lE - lLen + 1

  wdata <- data.frame(x=rep(NA,ndata),y=rep(NA,ndata))
  if(is.element("NAME",names(PW@data))) {
    cnames <- as.vector(PW@data$NAME)
  } else if(is.element("name",names(PW@data))) {
    cnames <- as.vector(PW@data$name)
  } else stop("Can't find name field.")

  gon.names <- character(totgon)
  k <- 1
  for(i in 1:ncountries){
    for (j in 1:ngon[i]){
      if(ngon[i]==1) gon.names[k] <- cnames[i]
      else gon.names[k] <- paste(cnames[i],":",j,sep="")
      cc <- data.frame(PW@lines[[i]]@Lines[[j]]@coords)
#      cat(i,j,k,lB[k],lE[k],"\n")
      wdata$x[lB[k]:lE[k]] <- cc[,1]
      wdata$y[lB[k]:lE[k]] <- cc[,2]
      k <- k+1
    }
  }
#  wdata <- round(wdata*pi/180*10^6)
# create list of default polyon names from country names list
  gLen <- rep(1,totgon)
  gE <- 2*(1:totgon)
  gB <- gE - 1
  gdata <- rep(NA,2*totgon)
  print(totgon)
  gdata[seq(1,2*totgon-1,by=2)] <- 1:totgon
  result <- list(line=list(x=wdata$x,y=wdata$y,B=lB,E=lE,length=lLen),
                 gon=list(names=gon.names,B=gB,E=gE,length=gLen,data=gdata),
                 ngon=ngon,ncountries=ncountries,cnames=cnames)
  result
}


