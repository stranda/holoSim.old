
suppressMessages({
    source("setup.R")
    library(rlecuyer, quiet=T)
    library(pbdMPI, quiet=T)
})

starttime <- Sys.time()
tmpdir <- "/tmp"
inputdir <- getwd()
resultsdir <- paste(inputdir,"results",sep="/")

init()

reps <- (comm.size()) * 1

refuge.types <- list( #treats refugia has an index to this list
	c(1),             # 1 one in lower left corner
	c(5),             # 2 bottom middle
	c(1,10),          # 3 both lower corners
        c(2,4,6,7,9),     # 4 half the bottom row
        c(5,100)          # 5 lower center, upper right
) 

treats <- expand.grid(
    dens.scale=0.001,
    refuges = c(2,3,4,5),
    refsize = c(1000),
    glacfront = c(20),
    marginal.decrease=c(1),
    shortshape=1,
    shortscale=c(1,0.5,0.25),
    nloc=100,
    longmean=c(1.5,2,2.5,3),
    mix=c(0.15,0.1,0.05,0.01,0),
    cpmut=10e-7,
    nmut=10e-5,
    reps=1
    )



simrep <- function(jid,treats)
    {
        print(jid)
        print(date())
	print(treats[jid,])
	res =  NULL
	x=jid
        refugia=refuge.types[[treats[x,"refuges"]]]
		res <- tryCatch({
                        run.rep(refugia=1,
                                      dens.scale=treats[x,"dens.scale"],
#                                      k=rep(250,100),
                                      refsize=treats[x,"refsize"],
                                      glac.front=treats[x,"glacfront"],
                                      marginal.decrease=treats[x,"marginal.decrease"],
                                      shortshape=treats[x,"shortshape"],
                                      shortscale=treats[x,"shortscale"],
                                      longmean=treats[x,"longmean"],
                                      mix=treats[x,"mix"],
                                      nmut=treats[x,"nmut"],
                                      cpmut=treats[x,"cpmut"],
                                      plotland=F,
#					nloc=10
                                      nloc=treats[x,"nloc"]
                                      )
                        
				},
                        error=function(err){err},
                	warning=function(warn){warn})
           list(treats=treats[x,],results=res)
    }

ret <- task.pull(1:dim(treats)[1], simrep, treats=treats)

if(comm.rank() == 0){ #sort of like a master
    save(file=paste0(resultsdir,"/",runif(1),"-simrep.rda"),ret)
    print(Sys.time()-starttime)
}

finalize()
