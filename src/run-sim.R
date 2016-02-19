#
# 
#
library(parallel)
source("setup.R")
reps=16

system.time(allreps <- mclapply(1:reps,mc.cores=9,function(x) {run.rep()}))
