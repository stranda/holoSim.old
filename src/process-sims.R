#
# take simulation outputs, convert to usable formats and run analyses on them
#
# also file the analyses away for future use and finally tally the treatment combinations
# to identify parameter space to investigate
#
library(RSQLite)
library(dplyr)
#library(tidyr)
library(compiler)
library(parallel)
source("setup.R")
source("analysis_functions.R")
resdir <- "results"
analyzeddir <- "analyzed"
#dbfile <- paste0("./",analyzeddir,"/","summary_stats.sqlite")
#filelst <- list.files(path=analyzeddir,pattern="summary_stats.sqlite")
#if (length(filelst)>0) resDB <- src_sqlite(dbfile) else  resDB <- src_sqlite(dbfile,create=T)
basefiles <- list.files(path=resdir, pattern="*.rda")
resfiles <- paste0(resdir,"/",basefiles)
nm <- system("hostname",intern=T)
if (nm[1]!="orchis") {resfiles <- rev(resfiles); cores <- 1} else {cores <- 4}

enableJIT(3)
### big loop to process the results

nlv <- 1:10
nsamp <- 30
for (f in 1:length(resfiles))
{
    if (!file.exists(paste0(analyzeddir,"/",basefiles[f])))
    {
        load(resfiles[f])
        gc()
        success=T
                                        #    ret <- ret[1:4]
        retlst <- #lapply (1:length(ret),function(i) #ret is the return object read in
            mclapply (1:length(ret),mc.cores=cores,function(i) #ret is the return object read in
            {
                print(paste(basefiles[f],"rep=",i,"of",length(ret),date()))
                ##first convert the stored landscape to gi, etc
                samp <- land2export(ret[[i]]$results$sample$land,
                                    ns=nsamp, pvec=NULL,nlocvec=nlv) #pvec can be used for sampling pops
                treat <- ret[[i]]$treats
                ares <- unlist(tryCatch(analysis.func(samp),error=function(e){NULL}))
                if (!is.null(ares))
                {
                    results <- treat[rep(1,length(ares)),]
                    results$statistic <- names(ares)
                    results$value <- c(ares)
                    results$samp.per.pop <- nsamp
                    results$numloci <- length(nlv)
                    results$computation <- basefiles[f]
                                        #        success=dbWriteTable(resDB$con,value=results,name="analysisResults",append=T)
                                        #        success
#                    write.csv(file=paste0("csvres/",gsub(".rda",paste0(i,".csv"),basefiles[f])),
#                              row.names=F,results)
                    results
                } else NULL
            })
        retdf <- do.call(rbind,retlst)
        save(file=paste0(analyzeddir,"/",basefiles[f]),retdf)
        rm(ret);gc()
    }

                                        #    if (sum(!unlist(flags))==0) {
#        file.copy(from=resfiles[f],to=paste0("procresults/",basefiles[f]))
#        file.remove(resfiles[f])
#    } else {stop("something went wrong with writing to a table")}

}
