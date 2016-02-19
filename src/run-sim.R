#
# 
#
library(parallel)
library(compiler)
enableJIT(3) #byte compile for speed
source("setup.R")
refuge.types <- list(c(1),c(5),c(1,10)) #treats refugia has an index to this list
treats <- expand.grid(
    refuges = c(1,3),
    refsize = 1000,
    glacfront = 20,
    marginal.decrease=1,
    shortshape=1,
    shortscale=1,
    nloc=100,
    longmean=3,
    mix=c(0,0.5),
    reps=1:20
    )

treats <- treats[1]

allreps <- mclapply(1:dim(treats)[1],mc.cores=9,function(x)
                    {
                        refugia=refuge.types[[treats[x,"refuges"]]]
                        rp <- run.rep(refugia=1,
                                      refsize=treats[x,"refsize"],
                                      glac.front=treats[x,"glacfront"],
                                      marginal.decrease=treats[x,"marginal.decrease"],
                                      shortshape=treats[x,"shortshape"],
                                      shortscale=treats[x,"shortscale"],
                                      longmean=treats[x,"longmean"],
                                      mix=treats[x,"mix"],
                                      nloc=treats[x,"nloc"]
                                      )
                        list(treats=treats[x],results=rp)
                    })
