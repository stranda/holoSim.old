library(ggplot2)
library(tidyr)


repin <- function(file="allreps-test1.rda")
{
    load(file)
                                        #convert to a more usable format
    m1 <- as.data.frame(t(sapply(allreps,function(x){c(x[[1]],unlist(x[[2]]))})))
    m2 <- as.data.frame(do.call(cbind,lapply(1:dim(m1)[2],function(i) unlist(m1[,i]))))
    names(m2) <- make.names(names(m1))
    classes <- c("refuges","refsize","glacfront", "marginal.decrease","shortshape","shortscale","nloc", "longmean","mix","reps")
    longdat <- gather_(data.frame(m2),'statistic','value',gather_cols=names(m2)[!(names(m2)%in%classes)])
}

dat <- rbind(repin("allreps-test1.rda"),repin("allreps-test2.rda"))

###########now do some plotting
###This first one emphasizes the differences among numbers/locations of refuges
### the refuge numbers map this way: 1=1 refuge lower left; 3=2 refuges lower left and right
###
dat$refuges <- as.factor(dat$refuges)
p <- ggplot(dat,aes(x=refuges,y=value)) + geom_boxplot(aes(group=refuges)) 
p <- p + facet_wrap(~statistic+mix, ncol=4, scales = "free_y")
p

##
## this one emphaizes the effects of gobs of LDD (mix=0.5)
##
dat$mix <- as.factor(dat$mix)
p <- ggplot(dat,aes(x=mix,y=value)) + geom_boxplot(aes(group=mix)) 
p <- p + facet_wrap(~statistic+refuges, ncol=4, scales = "free_y")
p
