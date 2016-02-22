#
#
# the basic idea is to simulate colonization from a dispersal kernel
# 
# start with a grid landscape.  The 'population' is the center of each grid
#   each time click for each population, pull a random number from the distribution
#   
#

make.pops <- function(h=100,extent=c(x0=0,y0=0,x1=100,y1=100),refugia=c(1))
{
    if (floor(sqrt(h))!=sqrt(h))
    {
        stop("must be square habitat")
    }
    grd <- expand.grid(row=1:sqrt(h),col=1:sqrt(h))
    grd <- grd[order(grd$row,grd$col),]
    grd$pop <- 1:dim(grd)[1]
    colwid <- abs(extent["x0"]-extent["x1"])/sqrt(h)
    rowheight <- abs(extent["y0"]-extent["y1"])/sqrt(h)
    colleft=min(grd$col*colwid)-(colwid/2)
    rowbot=min(grd$row*rowheight)-(rowheight/2)
    collocs <- seq(colleft,colleft+(sqrt(h))*colwid,colwid)
    rowlocs <- seq(rowbot,rowbot+(sqrt(h))*rowheight,rowheight)
    pops <- data.frame(t(sapply(grd$pop,function(x)
    {
        c(pop=x,
          row=grd$row[x],
          col=grd$col[x],
          x0=collocs[grd$col[x]],
          x1=collocs[grd$col[x]+1],
          y0=rowlocs[grd$row[x]],
          y1=rowlocs[grd$row[x]+1],
          arrived=NA,
          source=NA,
          quality=1
          )
    })))
    pops[pops$pop%in%refugia,"arrived"] <- 0
    pops[pops$pop%in%refugia,"source"] <- 0
    pops
}


#give it an x and y and find the pop.  if off the grid return NA
find.coords <- function(pops,x,y)
{
    p <- pops[with(pops,(x0<=x)&(x1>x)&(y0<=y)&(y1>y)),"pop"]
#    print(paste(x,y,p))
    if (length(p)>0) p else NA
}

#returns n numbers from the mixed distribution
rmixedpdf <- function(n=1,shortscale=1,shortshape=1,longmean=10,longvar=10,mix=0)
{
    weib <- rweibull(n,shape=shortshape,scale=shortscale)
    nrm <- rnorm(n,mean=longmean,sd=sqrt(longvar))
    choice <- runif(n)
    ifelse(choice<mix,nrm,weib)
}
                     

find.source <- function(pops,gen=0,shortscale=1,shortshape=1,longmean=10,longvar=10,mix=0)
{
    thetas <- 2*pi*runif(sum(is.na(pops$arrived)))
    dists <- rmixedpdf(length(thetas),shortscale=shortscale,shortshape=shortshape,longmean=longmean,longvar=longvar,mix=mix)
    pops$searchx <- NA
    pops$searchy <- NA
    pops$searchx[is.na(pops$arrived)] <- with(pops[is.na(pops$arrived),],(x1+x0)/2 )+ cos(thetas)*dists
    pops$searchy[is.na(pops$arrived)] <- with(pops[is.na(pops$arrived),],(y1+y0)/2 )+ sin(thetas)*dists

    popsrc <- apply(pops[is.na(pops$source),],1,function(x){
        find.coords(pops,x["searchx"],x["searchy"])
    })
    
    newpops <- popsrc!=pops[is.na(pops$source),"pop"] #these should be distances that are not local
    newpops[is.na(newpops)] <- FALSE

    newpops <- newpops & (popsrc %in% pops[!is.na(pops$arrived),"pop"])

    if (sum(newpops)>0)
    {
  #      print("there are new pops")
        pops[is.na(pops$source),]$source <- ifelse(newpops,popsrc,NA)
        pops[is.na(pops$arrived)&(!is.na(pops$source)),"arrived"] <- gen
    }
    
    pops[,c("pop","row","col","x0","x1","y0","y1","arrived","source","quality")]
}
##
##
##
run.sim <- function(h=100,
                    extent=c(x0=0,y0=0,x1=100,y1=100),
                    refugia=c(1),
                    maxit=500,
                    shortscale=1,
                    shortshape=1,
                    longmean=10,
                    longvar=10,
                    mix=0.01
                   )
{
    pops <- make.pops(h=h,extent=extent,refugia=refugia)
    
    cnt <- 1
    while ((cnt<=maxit)&(sum(is.na(pops$source))>0))
    {
        pops <- find.source(pops,gen=cnt,shortscale=shortscale,shortshape=shortshape,longmean=longmean,longvar=longvar)
        cnt <- cnt+1
    }
    pops
}

