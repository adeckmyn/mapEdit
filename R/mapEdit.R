### mapEdit
###
### edit maps from the maps and mapdata packages
### Alex Deckmyn
################################################
# run R 
# > m1 <- read.map("mapdata/src/worldHires")
# > m2 <- merge.map(m1)
# This second step gives you a "clean" set of lines and polygons,
# but takes quite a long time to run!


### INPUT: read source data

read.map <- function(name,scale=10^8){
  gon <- read.gon(name)
  gon <- parse.gon(gon)
  line <- read.line(name,scale=scale)
  line <- parse.line(line)
  list(gon=gon,line=line)
}

read.gon <- function(infile){
  data <- scan(paste(infile,'.gon',sep=''),na.strings='EOR')
  names <- read.table(paste(infile,'.name',sep=''),sep='\t',quote="",
                   stringsAsFactors=FALSE,
                   col.names=c('name','index'))$name
  list(data=data,names=names)
}

parse.gon <- function(gon){
### index the begin and end points of every gon
  NG <- length(gon$data)
  ngons <- sum(is.na(gon$data))
  gon.E <- which(is.na(gon$data))-1
  if(length(gon.E) != ngons) stop(paste('NGON error',ngons,length(gon.E)))
  gon.B <- c(1,gon.E[-ngons]+2)
  gon.NL <- gon.E-gon.B+1

  list(ngons=ngons,length=gon.NL,B=gon.B,E=gon.E,data=gon$data,names=gon$names)
}

read.line <- function(infile,scale=10^6){
  data <- read.table(paste(infile,'.line',sep=''),
                   col.names=c('x','y'),na.strings='EOR',fill=TRUE)
  data <- skip.doubles(data)
### first element of every line identifies the polygons it separates
  NL <- length(data$x)
  nlines <- sum(is.na(data$x))
  gonloc <- c(1,(which(is.na(data$x))+1))[1:nlines]
  line.L <- data$x[gonloc]
  line.R <- data$y[gonloc]
#
  data <- data[-gonloc,]

  data <- round(data*scale)

  list(x=data$x,y=data$y,nlines=nlines,L=line.L,R=line.R)
}

skip.doubles <- function(data){
### get rid of some repeated points
### they are rare, but exist in worldHires
### probably as a result of rounding long ago when saving in ASCII
  NL <- length(data$x)

### NA signifies the end of a line
### the first couple after NA is the identification  
  realpoints <- logical(NL)
  realpoints[2:NL] <- (!is.na(data$x)[2:NL] & !is.na(data$x)[1:(NL-1)])
  realpoints[1] <- FALSE

  dx <- abs(diff(data$x))
  dy <- abs(diff(data$y))
  realdiff <- realpoints[1:(NL-1)] & !is.na(dx)
### the "doubles" are points that are equal to the next point
  doubles <- realdiff & dx==0 & dy==0 
  print(paste('Found',sum(doubles),'double points (including mini-islands?)'))
### but some of these may be mini-islands! i.e. island loops with just 1 repeated point
### As the first entry for every line is in fact an identification
### so dx[(ddd-1)] can only be NA if this ID tag is identical to the first point
### very unlikely indeed
### to keep the "islands" we also require that there is a next point beyond the repeated point
### in fact, I don't believe there are any such "mini-island" dots in worldHires, 
### but in "world", even the very first line is one!
  realdoubles <- logical(NL-1) 
  realdoubles[2:(NL-2)] <- doubles[2:(NL-2)] & (realdiff[3:(NL-1)] | realdiff[1:(NL-3)] )
### actually, I know that doubles[1] and doubles[NL-1] are FALSE
### because these are not real points !
  realdoubles[1] <- doubles[1] & realdiff[2]
  realdoubles[(NL-1)] <- doubles[(NL-1)] & realdiff[(NL-2)]
  
  if(sum(realdoubles)>0){
    ddd <- (1:(NL-1))[realdoubles]
    print(paste('Removing',length(ddd),'double points')) 
    data <- data[-ddd,]
  }
  else {
    print(paste('No double points removed.')) 
   }
  data
}


parse.line <- function(ww=world){
### go through the lines data
### and create a list of begin and end indices etc.
### this routine is also run after changes to the lines (e.g. merging)
  NL <- length(ww$x)
  nlines <- sum(is.na(ww$x))
### locate beginning and end of lines
  line.E <- (1:NL)[is.na(ww$x)]-1
  if(length(line.E) != nlines) stop('NLINES error')
  line.B <- c(1,line.E[-nlines]+2)
### number of points in a line
  line.NP <- line.E-line.B+1

  list(x=ww$x,y=ww$y,nlines=nlines,L=ww$L,R=ww$R,B=line.B,E=line.E,length=line.NP)
}

flip.line <- function(loc,ww=world,mirror=TRUE){
### reverse the direction of a line
### mirroring (=adapt left and right region) isn't crucial 
### because these L and R are often switched anyway
### BUT we now want to make the data completely consistent
  newline <- ww$line
  newgon <- ww$gon

  Bl <- ww$line$B[loc]
  El <- ww$line$E[loc]
  newline$x[Bl:El] <- ww$line$x[El:Bl]
  newline$y[Bl:El] <- ww$line$y[El:Bl]
  if(mirror){
    newline$L[loc] <- ww$line$R[loc]
    newline$R[loc] <- ww$line$L[loc]
  }

  newgon$data[ww$gon$data==loc] <- -loc
  newgon$data[ww$gon$data==-loc] <- loc
  list(gon=newgon,line=newline)
}

#########################3
### LINE MERGING
### line.valence
### line.partner
### merge.line
### extend.line
### map.merge

line.valence <- function(loc,ww=world){
### calculate the valence of the end points of the polyline
### i.e. the number of line segments that end there.
### if it is 2, we can e.g. merge the lines.
  swap <- (loc<0)
  if(swap) loc <- -loc

  if(loc>ww$line$nlines) c(0,0)
  else {
    nb <- sum(ww$line$x==ww$line$x[ww$line$B[loc]] &
         ww$line$y==ww$line$y[ww$line$B[loc]], na.rm=TRUE )  
  
    ne <- sum(ww$line$x==ww$line$x[ww$line$E[loc]] &
         ww$line$y==ww$line$y[ww$line$E[loc]], na.rm=TRUE )  
### for a closed loop, return valence "0"
    if(ne==2 & nb==2 & 
       ww$line$x[ww$line$E[loc]]== ww$line$x[ww$line$B[loc]] &
       ww$line$y[ww$line$E[loc]]== ww$line$y[ww$line$B[loc]] ) 
     {
       ne <- 0
       nb <- 0
     }
    if(swap) c(ne,nb)
    else c(nb,ne)
  }
}

line.partner <- function(loc,end,ww=world){
### if a point has valence 2, find the two segments that can be merged
### if valence is 1, this is probably an error
### if valence is >2, it is a vertex and can not be merged.
  if(line.valence(loc,ww)[end]!=2) stop('Not valence 2')
  if(end==1) {
    X <- ww$line$x[ww$line$B[loc]]
    Y <- ww$line$y[ww$line$B[loc]]
  } else {
    X <- ww$line$x[ww$line$E[loc]]
    Y <- ww$line$y[ww$line$E[loc]]
  }
#  print(paste('matching',X,Y))

  bx <- ww$line$x[ww$line$B]
  by <- ww$line$y[ww$line$B]
  ex <- ww$line$x[ww$line$E]
  ey <- ww$line$y[ww$line$E]

  bb <- (1:ww$line$nlines)[bx==X & by==Y]
  ee <- (1:ww$line$nlines)[ex==X & ey==Y]

  if(end==1) bb=bb[bb!=loc]
  else ee <- ee[ee!=loc]
  if (length(bb)+length(ee)!=1) stop('length error')
  if(end==1){
    if(length(ee)==0) -bb
    else ee
  } else {
    if(length(bb)==0) -ee
    else bb
  }
}

merge.line <- function(l1,l2,ww=world,quiet=FALSE){
### function that actually merges two line segments
### this also requires changing the polygon definitions!
  if(!quiet) print(paste('Merge lines',l1,l2))
  if(l1<0) {l1 <- -l1;ww <- flip.line(l1,ww)}
  if(l2<0) {l2 <- -l2;ww <- flip.line(l2,ww)}


  if(ww$line$x[ww$line$E[l1]]!=ww$line$x[ww$line$B[l2]] |
     ww$line$y[ww$line$E[l1]]!=ww$line$y[ww$line$B[l2]] ) stop('No fit')

  if(l1==l2) stop('Trying to merge a line with itself? Maybe an island?')

  if(line.valence(l1,ww)[2]!=2 | line.valence(l2,ww)[1]!=2) stop('Wrong valence.')

  if(ww$line$L[l1]!=ww$line$L[l2] | ww$line$R[l1]!=ww$line$R[l2]) stop('different border regions!')

  nline <- ww$line
  ngon <- ww$gon
### FIX ME: if l1 or l2 == nline, l1+1 is meaningless !
  if(l2==l1+1) {
    nline$x <- ww$line$x[-(ww$line$E[l1]+1:2)]
    nline$y <- ww$line$y[-(ww$line$E[l1]+1:2)]
  }  else 
  if(l1==l2+1) {
      nline$x <- ww$line$x[c(1:ww$line$E[(l2-1)],NA,
                          ww$line$B[l1]:ww$line$E[l1],
                          (ww$line$B[l2]+1):ww$line$E[l2],NA,
                         ww$line$B[(l1+1)]:length(ww$line$x)  )]
      nline$y <- ww$line$y[c(1:ww$line$E[(l2-1)],NA,
                          ww$line$B[l1]:ww$line$E[l1],
                          (ww$line$B[l2]+1):ww$line$E[l2],NA,
                         ww$line$B[(l1+1)]:length(ww$line$x)  )]   
  }  else 
  if(l1<l2){
      nline$x <- ww$line$x[c(1:ww$line$E[l1],
                         (ww$line$B[l2]+1):ww$line$E[l2],
                         (ww$line$E[l1]+1):(ww$line$B[l2]-1),
                         (ww$line$E[l2]+2):length(ww$line$x))]
      nline$y <- ww$line$y[c(1:ww$line$E[l1],
                         (ww$line$B[l2]+1):ww$line$E[l2],
                         (ww$line$E[l1]+1):(ww$line$B[l2]-1),
                         (ww$line$E[l2]+2):length(ww$line$x))]
    } else {
       nline$x <- ww$line$x[c(1:ww$line$E[(l2-1)],NA,
                          ww$line$B[l1]:ww$line$E[l1],
                          (ww$line$B[l2]+1):ww$line$E[l2],NA,
                         (ww$line$B[(l2+1)]):(ww$line$E[(l1-1)]),NA,
                         ww$line$B[(l1+1)]:length(ww$line$x)  )]
       nline$y <- ww$line$y[c(1:ww$line$E[(l2-1)],NA,
                          ww$line$B[l1]:ww$line$E[l1],
                          (ww$line$B[l2]+1):ww$line$E[l2],NA,
                         (ww$line$B[(l2+1)]):(ww$line$E[(l1-1)]),NA,
                         ww$line$B[(l1+1)]:length(ww$line$x)  )]
    }

## the largest line number is now obsolete: renumber all following lines
  LL <- max(l1,l2)
  nline$L <- ww$line$L[-LL]
  nline$R <- ww$line$R[-LL]
  nline <- parse.line(nline)
### adapt the gon data !
  ngon$data <- ww$gon$data[abs(ww$gon$data)!=LL]
  xx <- (!is.na(ngon$data) & ngon$data>LL)
  ngon$data[xx] <- ngon$data[xx]-1
  xx <- (!is.na(ngon$data) & ngon$data< -LL)
  ngon$data[xx] <- ngon$data[xx]+1

  ngon <- parse.gon(list(data=ngon$data,names=ww$gon$names))
  list(gon=ngon,line=nline)
}

extend.line <- function(loc,ww=world){
  vv <- line.valence(loc,ww)
  if(vv[1]==2){
    nn <- line.partner(loc,1,ww)
    if(abs(nn)!=loc) ww <- merge.line(nn,loc,ww)                
  } else if(vv[2]==2){
    nn <- line.partner(loc,2,ww)
    if(abs(nn)!=loc) ww <- merge.line(loc,nn,ww)
  }
  ww
}

merge.map <- function(ww=world,nmax=Inf){
  count <- 0
  loc <- 1
  while(loc < ww$line$nlines) {
#    print(loc)
#    print(count)
    while(2 %in% line.valence(loc,ww) & count < nmax) {
      count <- count+1
      ww <- extend.line(loc,ww)
    }
    loc <- loc+1
  }
  print(paste('Merged',count,'polylines!'))
  ww
}

##################################3

get.line <- function(loc,ww=world){
  i1 <- ww$line$B[loc]
  i2 <- ww$line$E[loc]
  data.frame(x=ww$line$x[i1:i2],y=ww$line$y[i1:i2])
}

get.gon <- function(loc,ww=world){
  i1 <- ww$gon$B[loc]
  i2 <- ww$gon$E[loc]
  ww$gon$data[i1:i2]
}


replace.gon <- function(loc,newdata,newname,ww=world,corLR=FALSE){
  newgon <- ww$gon
  newgon$data <- c(ww$gon$data[1:ww$gon$E[(loc-1)]],NA,newdata,NA,
                            ww$gon$data[ww$gon$B[(loc+1)]:length(ww$gon$data)])
  newgon$names <- c(ww$gon$names[1:(loc-1)],newname,
                 ww$gon$names[(loc+1):ww$gon$ngons])

### in '$line', the border regions may change!
### We can only correct it for the new gon definition, assuming anti-clockwise orientation
### check this manually for all the line segments that may be affected
### i.e. those mentioned in newgon and oldgon
  newline <- ww$line
  if(corLR){
    for (ll in newdata){
      if (ll>0)  newline$L[ll] <- loc
      else newline$R[(-ll)] <- loc
    }
  }
  newgon <- parse.gon(newgon)
  list(gon=newgon,line=newline)
}



add.gon <- function(loc,newdata,newname,ww=world,corLR=FALSE){
  newgon <- ww$gon
  newgon$data <- c(ww$gon$data[1:ww$gon$E[(loc-1)]],NA,newdata,NA,
                            ww$gon$data[ww$gon$B[loc]:length(ww$gon$data)])
  newgon$names <- c(ww$gon$names[1:(loc-1)],newname,
                 ww$gon$names[loc:ww$gon$ngons])
  newgon$ngons <- ww$gon$ngons+1

### in 'line', the border regions move!

  newline <- ww$line
  if(corLR){
#    newline$L[newline$L>=loc] <-  newline$L[newline$L>=loc]+1
#    newline$R[newline$R>=loc] <-  newline$R[newline$R>=loc]+1
    newline$L <- newline$L + (newline$L>=loc)
    newline$R <- newline$R + (newline$R>=loc)
    for (ll in newdata){
      if (ll>0)  newline$L[ll] <- loc
      else newline$R[(-ll)] <- loc
    }
  }
  newgon <- parse.gon(newgon)
  list(gon=newgon,line=newline)
}

add.line <- function(loc,newx,newy,newL,newR,ww=world){
  newline <- ww$line
  newline$x <- c(ww$line$x[1:ww$line$E[(loc-1)]],NA,newx,NA,
              ww$line$x[ww$line$B[loc]:length(ww$line$x)])
  newline$y <- c(ww$line$y[1:ww$line$E[(loc-1)]],NA,newy,NA,
              ww$line$y[ww$line$B[loc]:length(ww$line$y)])
#  newline$nlines <- newline$nlines+1

  newline$L <- c(ww$line$L[1:(loc-1)],newL,ww$line$L[loc:ww$line$nlines])
  newline$R <- c(ww$line$R[1:(loc-1)],newR,ww$line$R[loc:ww$line$nlines])

  newline <- parse.line(newline)
# adapt the line numbers in the gon definitions (all lines > pos are moved up)
  newgon <- ww$gon
  newgon$data[!is.na(newgon$data) & newgon$data>=loc] <- newgon$data[!is.na(newgon$data) & newgon$data>=loc]+1
  newgon$data[!is.na(newgon$data) & newgon$data<= -loc] <- newgon$data[!is.na(newgon$data) & newgon$data<= -loc]-1
#  newgon <- parse.gon(newgon)
  list(gon=newgon,line=newline)
}

line.replace <- function(xy,loc,ww=world){
  NL <- ww$line$nlines
  newline <- ww$line
  if(loc==1){
    newline$x <-  c(xy$x,NA,
                  ww$line$x[ww$line$B[(loc+1)]:length(ww$line$x)])
    newline$y <-  c(xy$y,NA,
                  ww$line$y[ww$line$B[(loc+1)]:length(ww$line$y)])
  }
  else if (loc==NL) { 
    newline$x <-  c(ww$line$x[1:ww$line$E[(loc-1)]],NA,xy$x,NA)
    newline$y <-  c(ww$line$y[1:ww$line$E[(loc-1)]],NA,xy$y,NA)
  }
  else { 
    newline$x <-  c(ww$line$x[1:ww$line$E[(loc-1)]],NA,xy$x,NA,
                  ww$line$x[ww$line$B[(loc+1)]:length(ww$line$x)])
    newline$y <-  c(ww$line$y[1:ww$line$E[(loc-1)]],NA,xy$y,NA,
                  ww$line$y[ww$line$B[(loc+1)]:length(ww$line$y)])
  }
  newline <- parse.line(newline)
  ww$line <- newline
  ww
}


split.line <- function(pos,ww=world){
  if (pos %in% ww$line$B | pos %in% ww$line$E) stop('cannot split at vertex')
  lloc <- find.loc(pos,ww)
  newline <- ww$line
  newgon <- ww$gon

  newline$x <- c(ww$line$x[1:pos],NA,ww$line$x[pos:length(ww$line$x)])
  newline$y <- c(ww$line$y[1:pos],NA,ww$line$y[pos:length(ww$line$y)])
  newline$L <- c(ww$line$L[1:lloc],ww$line$L[lloc:ww$line$nlines])
  newline$R <- c(ww$line$R[1:lloc],ww$line$R[lloc:ww$line$nlines])
  newline <- parse.line(newline)

### the line may belong to 1 or 2 polygons (depends on L and R)
### but I know these are not always consistent
### so only use the sign of the line in the gon list
  newgon$data[!is.na(newgon$data) & newgon$data>lloc] <- 
             newgon$data[!is.na(newgon$data) & newgon$data>lloc]+1
  newgon$data[!is.na(newgon$data) & newgon$data< -lloc] <- 
             newgon$data[!is.na(newgon$data) & newgon$data< -lloc]-1

  lg <- ww$line$L[lloc]
  rg <- ww$line$R[lloc]

  noccur <- sum(abs(ww$gon$data)==lloc,na.rm=T)

  if(noccur>2) stop('the point is in more than 2 polygons')
  if(noccur<1) stop('the point does not appear in any polygon')
  
  for(i in 1:noccur){
### if there are 2 occurences -> the second time i -> i+1 !!!
### therefore, we recalculate the location
    gloc <- which(abs(newgon$data)==lloc)[i]
#(1:length(newgon$data))[!is.na(newgon$data) & lloc==abs(newgon$data)][i]
    if (newgon$data[gloc]>0) newgon$data <- c(newgon$data[1:gloc],lloc+1,
                                   newgon$data[(gloc+1):length(newgon$data)])
    else newgon$data <- c(newgon$data[1:(gloc-1)],-lloc-1,newgon$data[gloc:length(newgon$data)])
  }

  newgon <- parse.gon(newgon)
  list(gon=newgon,line=newline)
}

find.closest <- function(xy,ww=world){
  ddist <- abs(ww$line$x - xy$x)+abs(ww$line$y-xy$y)
#print(length(ddist))
  mindist <- min(ddist,na.rm=TRUE)
  print(paste('minimal distance:',mindist))
  print(paste('multiplicity:',sum(ddist==mindist,na.rm=TRUE)))
#  print(paste( (1:length(ww$line$x))[!is.na(dist)&dist==mindist])) 
  order(ddist)[1:10]

}

# which polyline does a particular point belong to...
find.loc <- function(i,ww=world){
  which(ww$line$B<=i & ww$line$E>=i)
}

# which polygons does a line belong to:
find.line <-  function(l,ww=world){
  zzz <- which(ww$gon$data == l)
  if(length(zzz)==0) return(NA)
  result <- zzz
  for(i in length(zzz)) result[i] <- which(ww$gon$B<=zzz[i] & ww$gon$E>=zzz[i])
  result
}



#######################
###########

### get a line out of another polyline set
get.poly <- function(n,dataset=world.unep){
    i1 <- (1:length(dataset$x))[is.na(dataset$x)][(n-1)]+1
    i2 <- (1:length(dataset$x))[is.na(dataset$x)][n]-1
#print(c(i1,i2))
    data.frame(x=dataset$x[i1:i2],y=dataset$y[i1:i2])
}

get.fullgon <- function(loc,ww=world){
  segm <- get.gon(loc,ww)
  ll <- get.line(abs(segm[1]),ww)
  if (segm[1]<0) ll <- ll[dim(ll)[1]:1,]
  contour <- ll

  for(i in segm[-1]){
    ll <- get.line(abs(i),ww)
    if (i<0) ll <- ll[dim(ll)[1]:1,]
    contour <- rbind(contour,ll[-1,])
  }
  contour
}

##############################################

line.borders <- function(loc,ww=world){
  for(i in loc) print(c(i,ww$line$L[abs(i)],ww$line$R[abs(i)]))
}
### return the line closest to given lat/lon
### default is that you click on the map to select the line!
closest.line <- function(data,xy=locator(1),quiet=FALSE){
  ll <- length(data$x)
  ddis <- abs(data$x-xy$x)+abs(data$y-xy$y)
  loc <- order(ddis)[1]
  if(!quiet)cat("Location: ",loc,"values:",data[loc,])
  i1 <- loc
  i2 <- loc
  while(!is.na(data$x[i1])) i1 <- i1-1
  while(!is.na(data$x[i2])) i2 <- i2+1
  data[(i1+1):(i2-1),]
    
}

get.map <- function(n,data=world.unep){
  xy <- locator(1)
  ll <- closest.line(xy,data)
  lines(ll,col=4)
  nmap <- ll
  for(i in 2:n){
    xy <- locator(1)
    ll <- closest.line(xy,data)
    lines(ll,col=4)
    nmap <- rbind(nmap,c(NA,NA),ll)
  }
  nmap
}
ap <- function(xy,w=w2,nx=1000,ny=nx,...){
  plot(w$line,type="o",xlim=c(-nx,nx)+xy$x,ylim=c(-ny,ny)+xy$y,...)
}

### CONSISTENCY CHECKS  

mirror.line <- function(loc,ww=world){
  i1 <- ww$line$L[loc]
  ww$line$L[loc] <- ww$line$R[loc]
  ww$line$R[loc] <- i1
  list(gon=ww$gon,line=ww$line)
}

check.gon <- function(gloc,ww=world){
### All countries should be on the "left" side of their borders (anti-clockwise)
### very useful to see if a gone is correctly formed
### notice that the original world database has many errors in L/R labels
### But that doesn't do much harm, apparantly.
  glines <- get.gon(gloc,data)
  for(lloc in glines){
    ll <- ww$line$L[abs(lloc)]
    rr <- ww$line$R[abs(lloc)]
    bx <- ww$line$x[ww$line$B[abs(lloc)]]
    by <- ww$line$y[ww$line$B[abs(lloc)]]
    ex <- ww$line$x[ww$line$E[abs(lloc)]]
    ey <- ww$line$y[ww$line$E[abs(lloc)]]
    if(lloc>0){
      print(c(ll,bx,lloc,by,rr))
      print(c(ll,ex,lloc,ey,rr))
    } else{
      print(c(rr,ex,lloc,ey,ll))
      print(c(rr,bx,lloc,by,ll))
    }
  }
}  

check.gon2 <- function(gloc,ww=world,correct=FALSE){
### the interior should always be on the left.
### also, all the line elements should fit together
### worldHires has a few tiny errors (UK, China)!
### This also finds even the smallest numerical inconsistencies
  glines <- get.gon(gloc,ww)
  N <- length(glines)
  result <- data.frame(matrix(NA,nrow=N,ncol=7))
  names(result) <- c("line","L","R","bx","by","ex","ey")
  if (N==1) return(if(correct) ww else NULL)
  for(i in 1:N){
    lloc <- glines[i]
    result$line[i] <- lloc
    if(lloc > 0){
      result$L[i] <- ww$line$L[abs(lloc)]
      result$R[i] <- ww$line$R[abs(lloc)]
      result$bx[i]  <- ww$line$x[ww$line$B[abs(lloc)]]
      result$by[i]  <- ww$line$y[ww$line$B[abs(lloc)]]
      result$ex[i]  <- ww$line$x[ww$line$E[abs(lloc)]]
      result$ey[i]  <- ww$line$y[ww$line$E[abs(lloc)]]
    }
    else {
      result$R[i] <- ww$line$L[abs(lloc)]
      result$L[i] <- ww$line$R[abs(lloc)]
      result$ex[i]  <- ww$line$x[ww$line$B[abs(lloc)]]
      result$ey[i]  <- ww$line$y[ww$line$B[abs(lloc)]]
      result$bx[i]  <- ww$line$x[ww$line$E[abs(lloc)]]
      result$by[i]  <- ww$line$y[ww$line$E[abs(lloc)]]
    }
  }
### check closure (usually OK)
  if(result$bx[1] != result$ex[N] | 
     result$by[1] != result$ey[N]) warning(paste("GON",gloc,": Bad closure 1"))
  for(i in 2:N) if(result$bx[i] != result$ex[(i-1)] | 
                   result$by[i] != result$ey[(i-1)]  )  warning(paste("GON",gloc,": Bad closure",i))
### check L/R placement:
  
  for(i in 1:N) if(result$L[i] != gloc) {
    warning(paste("GON",gloc,"line",i,"incorrectly oriented?"))
    if(correct) ww <- mirror.line(abs(glines[i]),ww)
  }
###
  if(correct) ww else NULL
}

###
# remove the ith point of a line
point.remove <- function(loc,line,ww=world){
  pl <- ww$line$B[line] + loc -1 
  ww$line$x <- ww$line$x[-pl]
  ww$line$y <- ww$line$y[-pl]
  ww$line$length[line] <- ww$line$length[line] - 1
  ww$line$B[(line+1):ww$line$nlines] <- ww$line$B[line:ww$line$nlines]-1
  ww$line$E[line:ww$line$nlines] <- ww$line$E[line:ww$line$nlines]-1
  ww
}

point.add <- function(xy,loc,line,ww=world){
  pl <- ww$line$B[line] + loc -1
  LL <- length(ww$line$x)
  ww$line$x <- c(ww$line$x[1:pl],xy[1],ww$line$x[(pl+1):LL])
  ww$line$y <- c(ww$line$y[1:pl],xy[2],ww$line$y[(pl+1):LL])
  ww$line$length[line] <- ww$line$length[line] + 1
  ww$line$B[(line+1):ww$line$nlines] <- ww$line$B[line:ww$line$nlines]+1
  ww$line$E[line:ww$line$nlines] <- ww$line$E[line:ww$line$nlines]+1
  ww
}

point.replace <- function(xy,loc,line,ww=world){
  pl <- ww$line$B[line] + loc -1
  ww$line$x[pl] <- xy[1]
  ww$line$y[pl] <- xy[2]
  ww
}

### find central points of all polygons
### this helps to find corresponding polygons in different sets
gonStat <- function(world){
  ngon <- world$gon$ngon
  result <- data.frame(x=rep(NA,ngon),y=rep(NA,ngon),L=rep(NA,ngon))
  for(i in 1:ngon){
    xy <- get.fullgon(i,world)
    result$x[i] <- mean(xy$x)
    result$y[i] <- mean(xy$y)
    result$L[i] <- dim(xy)[1]
  }
  result
}

gonFind<- function(gon,w1,w2){
  order( sqrt( (w1$x[gon]-w2$x)^2+(w1$y[gon]-w2$y)^2),na.last=T)[1:10]
}
getnames <- function(w1=world,x1=NULL,w2=wh,x2=NULL,begin=1,end=NULL){
  if(is.null(x1)) x1 <- gonStat(w1)
  if(is.null(x2)) x2 <- gonStat(w2)
  ngon2 <- w2$gon$ngon
  if(is.null(end)) end=ngon2
  for(i in begin:end){
    j=gonFind(i,x2,x1)[1]
    if(w1$gon$names[j] != w2$gon$names[i] ){
      xy2=get.fullgon(i,w2)
      xl=range(xy2$x)
      xlim=xl+diff(xl)*c(-1,1)*2
      yl=range(xy2$y)
      ylim=yl+diff(yl)*c(-1,1)*2

      plot(get.fullgon(i,w2),type="l",col=3,xlim=xlim,ylim=ylim)
      lines(get.fullgon(j,w1),type="l",col=2,lty=3,lwd=2)
      cat(j,w1$gon$names[j],i,w2$gon$names[i],"\n")
      x=readline()
      if(x=="s") return(w1)
      if(x!=""){
      print(x)
        if(x=="y") w1$gon$names[j]=w2$gon$names[i]
        else w1$gon$names[j]=x
      }
    }
  }
  w1
}


getnames2 <- function(w1=world,x1=NULL,w2=wh,x2=NULL,begin=1,end=NULL,swap=TRUE){
  if(is.null(x1)) x1 <- gonStat(w1)
  if(is.null(x2)) x2 <- gonStat(w2)
  ngon2 <- w2$gon$ngon
  if(is.null(end)) end=ngon2
  for(i in begin:end){
    j=gonFind(i,x2,x1)[1]
    if(w1$gon$names[j] != w2$gon$names[i] ){
      xy2=get.fullgon(i,w2)
      xl=range(xy2$x)
      xlim=xl+diff(xl)*c(-1,1)*2
      yl=range(xy2$y)
      ylim=yl+diff(yl)*c(-1,1)*2

      plot(get.fullgon(i,w2),type="l",col=3,xlim=xlim,ylim=ylim)
      lines(get.fullgon(j,w1),type="l",col=2,lty=3,lwd=2)
      cat(j,w1$gon$names[j],i,w2$gon$names[i],"\n")
      x=readline()
      if(x=="s") return(if(swap) w2 else w1)
      if(x!=""){
      print(x)
        if(x=="y") {if (!swap) w1$gon$names[j]=w2$gon$names[i] else w2$gon$names[i]=w1$gon$names[j] }
        else {if (!swap) w1$gon$names[j]=x else w2$gon$names[i]=x}
      }
    }
  }
  if(!swap) w1 else w2
}


#################3

# create (interpolated) points exactly on the central meridian
# PURPOSE: prepare for creation of a pacific-centered world map
meridian<-function(ww=world){
  i<-1
  k<-0
  while(i<length(ww$line$x)){
    if(!is.na(ww$line$x[i]) & !is.na(ww$line$x[(i+1)])){
      if(sign(ww$line$x[i])!=sign(ww$line$x[(i+1)])){
        xi <-0 
        yi <- ww$line$y[i] + (ww$line$y[(i+1)]-ww$line$y[i])/(ww$line$x[(i+1)]-ww$line$x[i]) * (xi-ww$line$x[i])
        cat(i,ww$line$x[i],xi,ww$line$x[i+1],"\n")
        cat(i,ww$line$y[i],yi,ww$line$y[i+1],"\n")
        ww <- point.add2(c(xi,yi),i+1,ww=ww)
        i <- i+1
        k<-k+1
      }
    }
    i <- i+1
  }
  cat("Added",k,"points\n")
  ww
}

# split lines that cross the central meridian (after inserting a split point at x=0)
# ATTENTION: after this, don't run map.merge !!! 
meridian.split<-function(ww=world){
  i<-2
  k<-0
  while(i<length(ww$line$x)){
    if(!is.na(ww$line$x[i]) & !is.na(ww$line$x[(i+1)]) & !is.na(ww$line$x[(i-1)])){
      if(ww$line$x[i]==0){
        cat("split at",i,ww$line$x[i],ww$line$x[i+1],"\n")
        ww <- split.line(i,ww)
        k<-k+1
      }
    }
    i <- i+1
  }
  cat("Splitted",k,"points\n")
  ww
}


point.add2 <- function(xy,pp,ww=world){
  pl <- pp-1
  loc <- which(ww$line$B<=pp & ww$line$E>=pp)
#  pl <- ww$line$B[line] + loc -1
  LL <- length(ww$line$x)
  ww$line$x <- c(ww$line$x[1:pl],xy[1],ww$line$x[(pl+1):LL])
  ww$line$y <- c(ww$line$y[1:pl],xy[2],ww$line$y[(pl+1):LL])
  ww$line$length[loc] <- ww$line$length[loc] + 1
  ww$line$B[(loc+1):ww$line$nlines] <- ww$line$B[(loc+1):ww$line$nlines]+1
  ww$line$E[loc:ww$line$nlines] <- ww$line$E[loc:ww$line$nlines]+1
  ww
}



### OUTPUT: routines to create the new source files 

export.map <- function(ww=world,outfile='world',scale=1,ndec=7){
# line data
  lfile <- paste(outfile,'.line',sep='')
  lsfile <- paste(outfile,'.linestats',sep='')

  lx <- round(ww$line$x * scale,ndec)
  ly <- round(ww$line$y * scale,ndec)
  system(paste('rm -f',lfile,lsfile))
  for(loc in 1:ww$line$nlines){
    write(paste(ww$line$L[loc],ww$line$R[loc]),file=lfile,append=TRUE)
    write(rbind('',format(
            rbind(lx[ww$line$B[loc]:ww$line$E[loc]],
                  ly[ww$line$B[loc]:ww$line$E[loc]]),
            nsmall=ndec)),
          file=lfile,append=T,ncolumns=3)

    write('EOR', file=lfile,append=T)
  }
# linestat
  system(paste('rm -f',lsfile))
  write(paste(ww$line$nlines,max(ww$line$length)),file=lsfile,append=TRUE)

# gon & name
  ind <- 1:ww$gon$ngons
  gfile <- paste(outfile,'.gon',sep='')
  gsfile <- paste(outfile,'.gonstats',sep='')
  nfile <- paste(outfile,'.name',sep='')
  system(paste('rm -f',gfile,nfile,gsfile))
  for(loc in 1:ww$gon$ngons){
### for exact match: 1 blank before the numbers
    write(paste('',ww$gon$data[ww$gon$B[loc]:ww$gon$E[loc]]),
          file=gfile,append=T,ncolumns=1)
    write('EOR', file=gfile,append=T)
    write(paste(ww$gon$name[loc],loc,sep='\t'),file=nfile,append=T)
  }
  write(paste(ww$gon$ngons,max(ww$gon$length)),file=gsfile,append=TRUE)
  
}

export.name <- function(ww=world,outfile='world'){
  nfile <- paste(outfile,'.name',sep='')
  system(paste('rm -f',nfile))

  for(loc in 1:ww$gon$ngons) write(paste(ww$gon$name[loc],loc,sep='\t'),file=nfile,append=T)
}

########################################3
### little bug: this function doesn't fix the L/R polygon data in world$line
### so you need to run 
remove.gon <- function(loc,ww=world){
  x1=ww$gon$B[loc]
  x2=ww$gon$E[loc]+1
  xl=ww$gon$length[loc]+1
  xt=ww$gon$ngons-1
  
  ww$gon$ngons <- ww$gon$ngons - 1
  ww$gon$names <- ww$gon$names[-loc]
  ww$gon$data <- ww$gon$data[-(x1:x2)]
  ww$gon$length <- ww$gon$length[-loc]
  ww$gon$B <- ww$gon$B[-loc]
  ww$gon$B[loc:xt] <- ww$gon$B[loc:xt] - xl
  ww$gon$E <- ww$gon$B + ww$gon$length - 1

  ww
}

remove.line <- function(loc,ww=world,quiet=TRUE){
  loc <- abs(loc)
  if (!quiet) print(paste("removing",loc))
  if(any(nna(ww$gon$data)==loc) | any(nna(ww$gon$data)== -loc)) stop("Line still used") 
  sect <- ww$line$B[loc]:(ww$line$E[loc]+1)
  ww$line$x <- ww$line$x[-sect]
  ww$line$y <- ww$line$y[-sect]
  ss <- ww$line$length[loc]+1
  ww$line$B <- ww$line$B[-loc] -
                     c(rep(0,loc-1),rep(ss,world$line$nlines-loc))
  ww$line$nlines <- ww$line$nlines-1
  ww$line$length <- ww$line$length[-loc]
  ww$line$E <- ww$line$B + ww$line$length - 1
  
  ww$line$L <- ww$line$L[-loc]
  ww$line$R <- ww$line$R[-loc]

  wx <- which(ww$gon$data > loc)
  ww$gon$data[wx] <- ww$gon$data[wx]-1
  wx <- which(ww$gon$data < -loc)
  ww$gon$data[wx] <- ww$gon$data[wx] + 1

  ww
}

orphan.lines <- function(ww){
  i <- 1
  while (i <= ww$line$nline) {
    zzz <- which(abs(ww$gon$data)==i)
    if(length(zzz)==0) ww <- remove.line(i,ww,quiet=FALSE)
    else i <- i+1
  }
  ww
}

remove.point <- function(i,ww=world){
  ll <- find.loc(i,ww)
  nl <- ww$line$nlines
  ww$line$x <- ww$line$x[-i]
  ww$line$y <- ww$line$y[-i]
  ww$line$E[ll:nl] <-  ww$line$E[ll:nl] - 1
  ww$line$B[(ll+1):nl] <-  ww$line$B[(ll+1):nl] - 1
  ww$line$length[ll] <-  ww$line$length[ll] - 1
  ww
}

###############################


