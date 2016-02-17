#
#
# run a single rep of a simulation.
#  this includes make a landscape, simulate deep history with simcoal and
#  then simulate 500 gens with rmetasim
#
#

run.rep <- function(refugia=c(1,3),
                    refsize=c(1000,1000),
                    dens.scale=0.0005,
                    samp.per.pop=24,
                    gens.per.epoch=50,
                    epochs=10,
                    splittime=120000,
                    cpmut=10e-08,
                    cpseq=500,
                    nloc=10,
                    nmut=10e-4,
                    h=100,
                    k=rep(100,h),
                    e=rep(0,h),
                    shortshape=1,shortscale=1,longmean=10,mix=0,
                    executable="fsc251"
                    )
    {

        l <- recolonizeLandscape(refs=refugia,sizeref=refsize,dens.scale=dens.scale,
                                 cpmut=cpmut,nssr=nloc,ssrmut=nmut,
                                 shortshape=shortshape,
                                 shortscale=shortscale,
                                 longmean=longmean,
                                 mix=mix)
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
#################
#################

        for (g in 1:epochs)
            {
                print(g)
                l <- landscape.simulate(l,gens.per.epoch)
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
        list(nucgi=gi,cpgi=multidna2genind(aln,mlst=T),cpaln=aln,land=sampland)
    }
