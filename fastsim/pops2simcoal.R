#
# takes a pops object and a bunch of other parameters and creates a character
# vector that contains a simcoal input file
#
pops2simcoal <- function(pops,popsizes=1000,popsamps=100)
{
    if (length(popsizes)==1)
        popsizevec <- rep(popsizes,dim(pops)[1]) else popsizevec <- popsizes
    if (length(popsamps)==1)
        popsampvec <- rep(popsamps,dim(pops)[1]) else popsampvec <- popsamps
    
    ret <- NULL
    ret <- "//Number of population samples (demes)"
    ret <- c(ret,dim(pops)[1])
    ret <- c(ret,"//Population effective sizes (number of genes)")
    ret <- c(ret,popsizevec)
    ret <- c(ret,"//Population sample sizes (number of genes)")
    ret <- c(ret,popsampvec)
}
