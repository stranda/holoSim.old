#
# this function is intended to create a simcoal object and run it
# then extract the genotypes and use them to initialize the rmetasm landscape
#
simcoal.init <- function(splittime=120000,
                         internode=0,
                         refugia=1,
                         refsize=1000,
                         cpmut=10e-06,
                         cpseq=500,
                         nloc=10,
                         nmut=10e-4,
                         executable="fsc251",
                         repID=NULL)
    {
        
        cp.par.file <- ifelse(is.null(repID),"cpDNA.par",paste0("cpDNA",repID,".par"))
        #run simcoal first for organelle
        simcoal(
            splittime=splittime,  #time back to single population
            internode = internode, #time between branching in a stepwise invasion
            marker.type= "mtDNA",# "mtDNA" or "usat"
            p=refugia,		# Number of populations
            n.popsize=refsize,	# Population size as the haploid number of chrom.
            n.samples=refsize,	# Number of individuals sampled from each population
            g.rate=0,	# Population growth rate
            mig.matrix=matrix(0,nrow=refugia,ncol=refugia),	# Migration matrix
            mt.bp=cpseq,	# mtDNA number of base pairs
            usat.nloci=0,	# usat number of loci
            mt.mut.rate=cpmut,# mtDNA mutation rate (per locus per generation)
            usat.mut.rate=0,	# usat mutation rate (per locus per generation)
            par.file=cp.par.file,   # Name of simcoal .par
            executable = "fsc251"
        )
        #run simcoal for nuclear ssr markers
        nuc.par.file <- ifelse(is.null(repID),"nuc.par",paste0("nuc",repID,".par"))
        simcoal(
            splittime=splittime,  #time back to single population
            internode = internode, #time between branching in a stepwise invasion
            marker.type="usat",# "mtDNA" or "usat"
            p=refugia,		# Number of populations
            n.popsize=refsize,	# Population size as the haploid number of chrom.
            n.samples=refsize,	# Number of individuals sampled from each population
            g.rate=0,	# Population growth rate
            mig.matrix=matrix(0,nrow=refugia,ncol=refugia),	# Migration matrix
            mt.bp=0,	# mtDNA number of base pairs
            usat.nloci=nloc,	# usat number of loci
            mt.mut.rate=0,# mtDNA mutation rate (per locus per generation)
            usat.mut.rate=nmut,	# usat mutation rate (per locus per generation)
            par.file=nuc.par.file,   # Name of simcoal .par
            executable = "fsc251"
        )

        ####now create some landscapes from the simcoal objects
        rl <- recolonizeLandscape(refs=1:refugia,sizeref=rep(100,refugia))
        rl$intparam$habitats <- refugia
        rl$loci <- NULL
        cpFileBase=gsub(".par","",cp.par.file)
        cpFile <- paste0(cpFileBase,"/",cpFileBase,"_1_1.arp")

        nucFileBase=gsub(".par","",nuc.par.file)
        nucFile <- paste0(nucFileBase,"/",nucFileBase,"_1_1.arp")
        ####this next landscape is only good for grafting in the loci and individuals into an actual landscape
        ####do not try to simulate.  don't worry if you try it will break
        rlnew <- landscape.coalinput(rl,npp=refsize,arlseq=cpFile,
                                     arlms=nucFile,
                                     seqsitemut=cpmut/cpseq,
                                     msmut=nmut)
#        unlink(cp.par.file)
#        unlink(nuc.par.file)

        rlnew
    }
                         
