simcoal <-
    function(
        splittime=120000,  #time back to single population
        internode = 0, #time between branching in a stepwise invasion
        marker.type,# "mtDNA" or "usat"
        p,		# Number of populations
        n.popsize,	# Population size as the haploid number of chrom.
        n.samples,	# Number of individuals sampled from each population
        g.rate,	# Population growth rate
        mig.matrix,	# Migration matrix
        mt.bp,	# mtDNA number of base pairs
        usat.nloci,	# usat number of loci
        mt.mut.rate,# mtDNA mutation rate (per locus per generation)
        usat.mut.rate,	# usat mutation rate (per locus per generation)
        par.file,   # Name of simcoal .par
        executable = "fsc251"
        )

{
    if (FALSE) #set to true for testing purposes (or just run the code block
        {
            splittime= 120000  #time back to single population
            internode = 0 #time between branching in a stepwise invasion
            marker.type ="usat"# "mtDNA" or "usat"
            p=3		# Number of populations
            n.popsize=rep(100,3)	# Population size as the haploid number of chrom.
            n.samples=rep(100,3)	# Number of individuals sampled from each population
            g.rate=0	# Population growth rate
            mig.matrix=matrix(0,ncol=p,nrow=p)	# Migration matrix
            mt.bp=100	# mtDNA number of base pairs
            usat.nloci=10	# usat number of loci
            mt.mut.rate=10e-04  # mtDNA mutation rate (per locus per generation)
            usat.mut.rate=10e-04	# usat mutation rate (per locus per generation)
            par.file ="test.par"  # Name of simcoal .par
            executable = "fsc251"
        }
                                        #browser()
    
    if(!any(marker.type==c("mtDNA","usat")))
        stop("marker.type must be either \"mtDNA\" or \"usat\"")
    
    if((p!=length(n.popsize))|(p!=length(n.samples)))
        stop("Something is wrong with dimensions")
    
    
    cat("//Input parameters for the coalescence simulation program : simcoal\n",file=par.file)
    cat(p,"samples to simulate\n",file=par.file,append=T)
    cat("//Deme sizes (haploid number of genes)\n",file=par.file,append=T)
    cat(paste0(n.popsize,sep="\n",collapse=""),file=par.file,append=T)
    cat("//Sample sizes\n",file=par.file,append=T)
    cat(paste0(n.samples,sep="\n",collapse=""),file=par.file,append=T)
    cat("//Growth rates\n",file=par.file,append=T)
    cat(paste0(rep(0,p),sep="\n",collapse=""),file=par.file,append=T)
    cat("//Number of migration matrices : If 0 : No migration between demes\n",file=par.file,append=T)
    if (p==1)
        {
            cat(0,"\n",file=par.file,append=T)
        }
    else
        {
            if ((max(mig.matrix) >= 10e-5)&(internode==0))
                {
                    cat(1,"\n",file=par.file,append=T)
                    cat("//Migration rates matrix 0 : \n",file=par.file,append=T)
                    write.table(matrix(sprintf("%f",mig.matrix),nrow=nrow(mig.matrix)),file=par.file,append=T,quote=F,row=F,col=F)
                } else {
                    cat(0,"\n",file=par.file,append=T)
                }
        }
    cat("//Historical event: time, source, sink, proportion of migrants, new deme size, new growth rate, new migration matrix\n",file=par.file,append=T)
    
    hist.events <- p-1
    
    cat(paste(hist.events,"historical events\n"),file=par.file,append=T)
    
    if ((internode*hist.events)>splittime)
        {
            internode <- floor(splittime/hist.events)
            message("internode length X number of internodes greater than total time depth")
        }
    
    if (p>1)
        {
            cnt <- (hist.events-1)
            for (i in 1:(hist.events))
                {
                    tau_i <- splittime - floor(cnt * internode)
                    cat(paste(as.integer(tau_i),""),file=par.file,append=T)
                    cat(paste(i,"0 1 2 0 0\n"),file=par.file,append=T)
                cnt <- cnt-1
                }
        }
    
    cat("//Number of independent (unlinked) chromosomes, and \"chromosome structure\" flag:  0 for identical structure across chromosomes, and 1 for different structures on different chromosomes.\n",file=par.file,append=T)
    if(marker.type=="mtDNA")
        cat(c(1,0),"\n",file=par.file,append=T)
    else
        cat(c(usat.nloci,0),"\n",file=par.file,append=T)
    cat("//Number of contiguous linkage blocks in chromosome 1:\n1\n",file=par.file,append=T)
    cat("//Per Block: Data type, No. of loci, Recombination rate to the right-side locus, plus optional parameters ***see detailed explanation here***\n",file=par.file,append=T)
    if(marker.type=="mtDNA")
        cat("DNA",mt.bp,0,mt.mut.rate,0,"\n",file=par.file,append=T)
    else
        cat("MICROSAT",1,0.2,usat.mut.rate,0,0,"\n",file=par.file,append=T)
                                        #AES changed to simcoal2 3/30/09
###set up for unique reps
    tmpres=system(paste0(executable," -n 1 -S -i ",par.file),T)
}
