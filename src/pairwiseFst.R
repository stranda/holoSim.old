#function name is self-evident.  USes the non-bias corrected wright function
#
landscape.Fst.pairwise <- function(rland)
    {
        #find populations that exist
        poplst <- as.numeric(names(table(landscape.populations(rland))))
        #make matrix
        retlst <- vector("list",(rland$intparam$habitats^2-rland$intparam$habitats)/2)
        rlcnt <- 1
        df <- expand.grid(p1=1:rland$intparam$habitats,
                          p2=1:rland$intparam$habitats)
        df <- df[df[,1]!=df[,2],]
        df <- unique(t(apply(df,1,sort)))

        retdf <- apply(df,1,function(x)
                       {
                           row=x[1]
                           col=x[2]
                           if ((row %in% poplst)&(col%in%poplst))
                               {
                                   sl <- landscape.sample(rland,pvec=c(row,col))
                                   rowMeans(landscape.Fst(sl),na.rm=T)
                                } else {rep(NA,rland$locnum)}
                       })

    }



landscape.Fst.pairwise <- function(rland)
    {
        #find populations that exist
        poplst <- as.numeric(names(table(landscape.populations(rland))))
        #make matrix
        retlst <- vector("list",(rland$intparam$habitats^2-rland$intparam$habitats)/2)
        rlcnt <- 1
        for (row in 2:dim(retmat)[1])
            for (col in 1:(row-1))
                {
                    print(rlcnt)
#                    if ((row %in% poplst)&(col%in%poplst))
#                        {
#                            sl <- landscape.sample(rland,pvec=c(row,col))
#                            retlst[[rlcnt]] <- list(pop1=row,pop2=col,Fst.loci=rowMeans(landscape.Fst(sl),na.rm=T))
#                        }
                    rlcnt <- rlcnt+1
                }

    }
