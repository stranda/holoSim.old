library(ggplot2)
library(dplyr)
analyzedir <- "analyzed"
infiles <- list.files(path=analyzedir,pattern="*.rda")
fullnames <- paste0(analyzedir,"/",infiles)

dat <- do.call(rbind,lapply(fullnames,function(x) {load(x);print(x);print(dim(retdf));retdf}))

#########make a data frame that shows the numbers of reps per treatment combinations
#########creates a file called 'rep-summary.csv'
reps <- select(dat,dens.scale,refuges , refsize ,glacfront,marginal.decrease,shortshape,shortscale,nloc,longmean,mix,cpmut,nmut,samp.per.pop,numloci,computation,reps) %>% distinct() %>% select(-reps,-computation) %>% table() %>% as.data.frame()
write.csv(file="rep-summary.csv",row.names=F,reps)

###########now do some plotting

mixes <- unique(dat$mix)
refs <- unique(dat$refuges)
longmn <- unique(dat$longmean)
sscl <-  unique(dat$shortscale)
frnt <- unique(dat$glacfront)
mrgs <- unique(dat$marginal.decrease)


pdf(file="plot-mix.pdf",paper="special",width=10,height=16)
for (gl in frnt)
    for (mrg in mrgs)
for (rf in refs)
    for(lmn in longmn)
        for(ss in sscl)
        {
            p <-filter(dat,refuges==rf,longmean==longmn,shortscale==sscl) %>%
                ggplot(aes(x=mix,y=value))
            p <- p + facet_wrap(~statistic,scales="free",ncol=4)
            p <- p+geom_boxplot(aes(group=mix))
            p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
            p <- p+labs(title=paste0("Refuges: ",rf,", longmean: ",lmn,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
            print(p)
        }
dev.off()

pdf(file="plot-refuges.pdf",paper="special",width=10,height=16)

for (gl in frnt)
    for (mrg in mrgs)
for (mx in mixes)
    for(lmn in longmn)
        for(ss in sscl)
        {
            p <-filter(dat,mixes=mx,longmean==longmn,shortscale==sscl) %>%
                ggplot(aes(x=refuges,y=value))
            p <- p + facet_wrap(~statistic,scales="free",ncol=4)
            p <- p+geom_boxplot(aes(group=refuges))
            p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
            p <- p+labs(title=paste0("Mix: ",mx,", longmean: ",lmn,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
            print(p)
        }
dev.off()


pdf(file="plot-longmean.pdf",paper="special",width=10,height=16)
for (gl in frnt)
    for (mrg in mrgs)
        for (mx in mixes)
            for(rf in refs)
                for(ss in sscl)
                {
                    p <-filter(dat,mixes=mx,refuges==rf,shortscale==sscl) %>%
                        ggplot(aes(x=longmean,y=value))
                    p <- p + facet_wrap(~statistic,scales="free",ncol=4)
                    p <- p+geom_boxplot(aes(group=longmean))
                    p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
                    p <- p+labs(title=paste0("Mix: ",mx,", refuges: ",rf,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
                    print(p)
                }
dev.off()

