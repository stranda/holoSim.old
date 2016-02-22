#
# 
#
library(parallel)
library(compiler)
enableJIT(3) #byte compile for speed
source("setup.R")
refuge.types <- list(c(1),c(5),c(1,10)) #treats refugia has an index to this list
treats <- expand.grid(
    refuges = c(1),
    refsize = 1000,
    glacfront = 20,
    marginal.decrease=1,
    shortshape=1,
    shortscale=1,
    nloc=1,
    longmean=3,
    mix=c(0),
    reps=1
    )

#treats <- treats[1:4,]

id <- paste0(date(),'-',round(runif(1,0,10000000)))
repSamples <- mclapply(1:dim(treats)[1],mc.cores=8,function(x)
                    {
                        refugia=refuge.types[[treats[x,"refuges"]]]
                        rp <- run.rep(refugia=1,
#                                      k=rep(250,100),
                                      refsize=treats[x,"refsize"],
                                      glac.front=treats[x,"glacfront"],
                                      marginal.decrease=treats[x,"marginal.decrease"],
                                      shortshape=treats[x,"shortshape"],
                                      shortscale=treats[x,"shortscale"],
                                      longmean=treats[x,"longmean"],
                                      mix=treats[x,"mix"],
                                      cpmut=0,
                                      nmut=0,
                                      plotland=F,
                                      nloc=treats[x,"nloc"]
                                      )
                        list(treats=treats[x,],results=rp)
                    })
save("results/repSamples-",id,".rda",repSamples)
