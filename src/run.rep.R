###
###
### run a single rep of a simulation.
###  this includes make a landscape, simulate deep history with simcoal and
###  then simulate 500 gens with rmetasim
###
###

#################
getRepSample <- function(l,glac.front=20,samp.per.pop=100,plotland=F)
		{
	      if (plotland)
                    {
			require(lattice)
			crd=landscape.popcoord(l)
                        print(paste("glac front",glac.front))
                        tbl <- data.frame(table(landscape.populations(l)))
                        names(tbl)[1] <- "pop"
                        plotpops <- merge(tbl,crd)
                        plt <- levelplot(Freq~col+row,data=plotpops,xlim=c(0,max(crd$row)+1),
                                                ylim=c(0,max(crd$col)+1),
                                                panel=function(...){panel.levelplot(...)
                                                                    panel.abline(h=glac.front)
                                                                }
                                                )
                    }
		

	sampland <- landscape.sample(l,ns=samp.per.pop)

        gi <- landscape.make.genind(sampland) #ignores cpdna

        seqs <- landscape.locus.states(sampland,1) #start to get the cpdna seqs worked out
        seqind <- sapply(sampland$individuals[,7],function(ai){seqs$state[which(seqs$aindex==ai)]})
        aln <- new("multidna")
        aln@dna=list(cpDNA=as.DNAbin(do.call(rbind,strsplit(tolower(seqind),""))))
        aln@labels=as.character(sampland$individuals[,4])
        aln@ind.info=data.frame(sampland$individuals[,c(1,4,5,6)])
        names(aln@ind.info) <- c("state","ID","MatID","PatID")
        #return gi, gi for cp, alignmnent for cp and the sampled landscape
        if (plotland)
            simres <- list(nucgi=gi,cpgi=multidna2genind(aln,mlst=T),cpaln=NULL,land=sampland,plots=plots)
        else
            simres <- list(nucgi=gi,cpgi=multidna2genind(aln,mlst=T),cpaln=NULL,land=sampland)
#        save(file="current-simrep.rda",simres)
#        analysis.func(simres,doplots=F)
        simres
    }



run.rep <- function(refugia=c(1,2),
                    refsize=c(1000,1000),
                    dens.scale=0.001,
                    samp.per.pop=100,
                    gens.per.epoch=50,
                    epochs=10,
                    splittime=120000,
                    cpmut=10e-06,
                    cpseq=500,
                    nloc=10,
                    nmut=10e-4,
                    h=100,
                    k=rep(1000,h),
                    e=rep(0,h),
                    shortshape=1,
                    shortscale=1,
                    longmean=3,
                    mix=0,
                    executable="fsc251",
                    glac.retreat.per.epoch=1,
                    decay.dist=3,
                    glac.front = sqrt(h)+decay.dist, #default is no glacier effect
                    marginal.decrease=1,             #default is no edge effect
                    dumbinit=F,
                    plotland=F
                    )
    {
        if (FALSE)
            {
               refugia=c(1,2)
               refsize=c(1000,1000)
               dens.scale=0.0005
               samp.per.pop=24
               gens.per.epoch=50
               epochs=10
               splittime=120000
               cpmut=10e-08
               cpseq=500
               nloc=10
               nmut=10e-4
               h=100
               k=rep(250,h)
               e=rep(0,h)
               shortshape=1
               shortscale=1
               longmean=10
               mix=0
               executable="fsc251"
               glac.front =15
               glac.retreat.per.epoch=1
               decay.dist=3
               marginal.decrease=0.5
               plotland=T
            }
        l <- recolonizeLandscape(refs=refugia,sizeref=refsize,dens.scale=dens.scale,
                                 k=k,e=e,
                                 cpmut=cpmut,nssr=nloc,ssrmut=nmut,
                                 shortshape=shortshape,
                                 shortscale=shortscale,
                                 longmean=longmean,
                                 mix=mix)
        if (!dumbinit)
        {
            fscl <- simcoal.init(refugia=length(refugia),
                                 refsize=refsize, cpseq=cpseq, cpmut=cpmut, nloc=nloc, nmut=nmut,
                                 executable=executable,
                                 repID=paste0(round(runif(1)*100000),round(runif(1)*100000))
                                 )
#### #do simcoal simulation and return a bastardized object
################## meld the fastsimcola landscape with the working rmetasim landscape.  take care
################## to make sure that refugia work
##################
            l$intparam$locusnum <- length(fscl$loci)
            l$loci <- fscl$loci
            l$individuals <- fscl$individuals
            pops <- data.frame(oldpops=landscape.populations(l))
            pops$newpops <- NA
            for (p in 1:length(refugia))
            {
                pops$newpops[pops$oldpops==p] <- refugia[p]
            }
            l$individuals[,1] <- (pops$newpops-1)*2
            tmpseq <- l$loci[[1]]$alleles[[1]]$state
            len <- nchar(tmpseq)
            tmpseq <- gsub("N","",tmpseq)
            newlen <- nchar(tmpseq)
            addseq <- paste(sample(c("A","C","G","T"),len-newlen,replace=T),collapse='')
            l$loci[[1]]$alleles <- lapply(l$loci[[1]]$alleles,
                                          function(x)
                                          {
                                              x$state <-  gsub("N","",x$state)
                                              x$state <- paste0(x$state,addseq)
                                              x
                                          })
        }
#################
################# actually do the rmetasim simulations
#################


	resvec <- vector("list",epochs+1)
	resvec[[1]] <- getRepSample(l,glac.front=glac.front,plotland=plotland)
        print("runnning simulation")
        for (g in 1:epochs)
            {
                print(g)
                ##set carry to correspond to glaciated front
                crd <- landscape.popcoord(l)
                crd$carry.fact <- 1.0
#                crd$carry.fact[(glac.front-crd$row)>decay.dist] <- 1.0  
                crd$carry.fact[((glac.front)-crd$row)<=decay.dist] <- (1/decay.dist)*((glac.front)-crd$row[((glac.front)-crd$row)<=decay.dist])
                crd$carry.fact[crd$row>=(glac.front)] <- 0
                crd$carry.fact[crd$col%in%c(1,max(crd$col))] <- crd$carry.fact[crd$col%in%c(1,max(crd$col))]*marginal.decrease
                l$demography$epochs[[1]]$Carry <- k*crd$carry.fact
                print(unique(crd[,c("row","carry.fact")]))
                l <- landscape.simulate(l,gens.per.epoch)

             	resvec[[g+1]] <- getRepSample(l,glac.front=glac.front,plotland=plotland)   
                    
                glac.front <- glac.front+glac.retreat.per.epoch
            }

#################
#################
	resvec
}
