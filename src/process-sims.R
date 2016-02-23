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
dbfile <- paste0("./",analyzeddir,"/","summary_stats.sqlite")
filelst <- list.files(path=analyzeddir,pattern="summary_stats.sqlite")
if (length(filelst)>0) resDB <- src_sqlite(dbfile) else  resDB <- src_sqlite(dbfile,create=T)
basefiles <- list.files(path=resdir, pattern="*.rda")
resfiles <- paste0(resdir,"/",basefiles)

enableJIT(3)
### big loop to process the results
for (f in 1:length(resfiles))
{
    load(resfiles[f])
    success=T
    flags <- mclapply (1:length(ret),mc.cores=4,function(i) #ret is the return object read in
    {
        print(paste(basefiles[f],"rep=",i,"of",length(ret)))
        ##first convert the stored landscape to gi, etc
        samp <- land2export(ret[[i]]$results$sample$land)
        treat <- ret[[i]]$treats
        ares <- unlist(analysis.func(samp))
        results <- treat[rep(1,length(ares)),]
        results$statistic <- names(ares)
        results$value <- c(ares)
        results$computation <- basefiles[f]
        success=dbWriteTable(resDB$con,value=results,name="analysisResults",append=T)
        success
    })
    if (sum(!unlist(flags))==0) {
        file.copy(from=resfiles[f],to=paste0("procresults/",basefiles[f]))
        file.remove(resfiles[f])
    } else {stop("something went wrong with writing to a table")}
    rm(ret);gc()
}
