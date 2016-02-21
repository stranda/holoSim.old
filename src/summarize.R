library(ggplot2)
library(tidyr)


load("allreps-test1.rda")

#convert to a more usable format

m1 <- as.data.frame(t(sapply(allreps,function(x){c(x[[1]],unlist(x[[2]]))})))
m2 <- as.data.frame(do.call(cbind,lapply(1:dim(m1)[2],function(i) unlist(m1[,i]))))
names(m2) <- make.names(names(m1))
classes <- c("refuges","refsize","glacfront", "marginal.decrease","shortshape","shortscale","nloc", "longmean","mix","reps")
gather_(data.frame(m1),'statistic','value',gather_cols=names(m1)[!(names(m1)%in%classes)])

