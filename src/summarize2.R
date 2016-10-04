library(ggplot2)
library(dplyr)
analyzedir <- "analyzed"
infiles <- list.files(path=analyzedir,pattern="*.rda")
fullnames <- paste0(analyzedir,"/",infiles)

dat <- do.call(rbind,lapply(fullnames,function(x) {
    
    load(x)
#    print(x);print(dim(retdf))
##the next few lines are a hack to pull out analyses that failed    
    if (sum(is.na(as.numeric(as.character(retdf$mix))))>0)
    {
        retdf <- retdf[retdf$mix==FALSE,]
        print("there was a problem analyzing this simulation set")
        print(x);print(dim(retdf))
    }
    retdf
}))

print(paste(dim(dat)[1],"statistic/simulation_rep combinations"))


#########make a data frame that shows the numbers of reps per treatment combinations
#########creates a file called 'rep-summary.csv'
reps <- select(dat,dens.scale,refuges , refsize ,glacfront,marginal.decrease,shortshape,shortscale,nloc,longmean,mix,cpmut,nmut,samp.per.pop,numloci,computation,reps) %>% distinct() %>% select(-reps,-computation) %>% table() %>% as.data.frame()#%>%filter(Freq>0)
write.csv(file="rep-summary.csv",row.names=F,reps)

###########now do some plotting

#establish ranges of input parameters in dat
mixes <- unique(dat$mix)
refs <- unique(dat$refuges)
longmn <- unique(dat$longmean)
sscl <-  unique(dat$shortscale)
frnt <- unique(dat$glacfront)
mrgs <- unique(dat$marginal.decrease)

#subset ranges to focus on certian combinations
frnt <- max(frnt) #only take a large value for glacfront. Means the glaciers never came south; habitat always good
mrgs <- max(mrgs) #highest marginal.decrease (named backwards) will be equal to 1 and means no increase or decrease
#refs <- refs[refs!=3] #bottom left only and bottom and embedded
mixes <- mixes[mixes<=0.1]
sscl <- sscl[sscl==0.5]

###now subset the data to only include the ranges identified above
dat <- select(dat,-reps)%>%
    filter(refuges%in%refs,longmean%in%longmn,shortscale%in%sscl,
           glacfront%in%frnt,marginal.decrease%in%mrgs,mix%in%mixes)



################################################
################################################



stats <- unique(dat$statistic)
stats <- stats[-grep("pval",stats)]
tdf <- filter(dat,statistic%in%stats,glacfront==20,marginal.decrease==1,shortscale==0.5,longmean==1.5,refuges==1)
if (dim(tdf)[1]>0)
{
    p <- tdf %>% ggplot(aes(x=mix,y=value))
    p <- p+geom_boxplot(aes(group=mix))
    p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
    p <- p + facet_wrap(~statistic,ncol=3,scales="free_y")
    p <- p + labs(y="Summary statistic value",x="Proportion of dispersal events considered Long-Distance")
                                        #p <- p+labs(title=paste0("Statistic: ",stat))
    png("highlight-effects.png",width=800,height=1050)
    print(p)
    dev.off()

    png("summary-stat-placeholder.png",width=500,height=300)
    p = tdf %>% filter(statistic%in%c("ALLELE_STATS.slope","FIT_DISTRIBUTION_WEIBULL.scale"))%>% ggplot(aes(x=mix,y=value))
     p <- p+geom_boxplot(aes(group=mix))
    p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
    p <- p + facet_wrap(~statistic,ncol=3,scales="free_y")
    p <- p + labs(y="Summary statistic value",x="Proportion of dispersal events considered Long-Distance")
                                        #p <- p+labs(title=paste0("Statistic: ",stat))
    
    print(p)
    dev.off()
    
    pdf("highlight-effects.pdf",width=10,height=16)
    print(p)
    dev.off()
}





pdf(file="plot-mix.pdf",paper="special",width=10,height=16)
for (gl in frnt)
    for (mrg in mrgs)
        for (rf in refs)
            for(lmn in longmn)
                for(ss in sscl)
                {
                    p <-filter(dat,refuges==rf,glacfront==gl,marginal.decrease==mrg,longmean==lmn,shortscale==ss) %>%
                        ggplot(aes(x=mix,y=value))
                    p <- p + facet_wrap(~statistic,scales="free_y",ncol=4)
                    p <- p+geom_boxplot(aes(group=mix))
                    p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
                    p <- p+labs(title=paste0("Refuges: ",rf,", longmean: ",lmn,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
                    print(p)
                }
dev.off()


print("finished plotting mix")

pdf(file="plot-refuges.pdf",paper="special",width=10,height=16)

for (gl in frnt)
    for (mrg in mrgs)
        for (mx in mixes)
            for(lmn in longmn)
                for(ss in sscl)
                {
                    print(c(gl,mrg,mx,lmn,ss))
                    tdf <- filter(dat,mix=mx,longmean==lmn,shortscale==ss,
                                  marginal.decrease==mrg,glacfront==gl)
                    if (dim(tdf)[1]>0)
                    {
                        tdf$refuges <- as.factor(tdf$refuges)
                        p <- tdf  %>% ggplot(aes(x=refuges,y=value))
                        p <- p + facet_wrap(~statistic,scales="free_y",ncol=4)
                        p <- p+geom_boxplot(aes(group=refuges))
                        p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
                        p <- p+labs(title=paste0("Mix: ",mx,", longmean: ",lmn,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
                        print(p)
                    } else {print("combo with no reps")}
                }
dev.off()

print("finished plotting refuges")

pdf(file="plot-longmean.pdf",paper="special",width=10,height=16)
for (gl in frnt)
    for (mrg in mrgs)
        for (mx in mixes)
            for(rf in refs)
                for(ss in sscl)
                {
                    tdf <- filter(dat,mix=mx,refuges==rf,shortscale==ss,glacfront==gl,marginal.decrease==mrg)
                    if (dim(tdf)[1]>0)
                    {
                        p <- tdf %>%
                            ggplot(aes(x=longmean,y=value))
                        p <- p + facet_wrap(~statistic,scales="free_y",ncol=4)
                        p <- p+geom_boxplot(aes(group=longmean))
                        p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
                        p <- p+labs(title=paste0("Mix: ",mx,", refuges: ",rf,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
                        print(p)
                    }
                }
dev.off()


pdf("mix-longmean.pdf")
for (gl in frnt)
    for (mrg in mrgs)
        for(rf in refs)
            for(ss in sscl)
            {
                tdf <- filter(dat,refuges==rf,shortscale==ss,glacfront==gl,marginal.decrease==mrg) %>%
                    group_by(longmean,mix,statistic) %>% summarise(value=mean(value))
                if (dim(tdf)[1]>0)
                {
                    p <- tdf %>% ggplot(aes(x=longmean,y=mix,z=value))
                    p <- p+geom_contour()
                    p <- p + facet_wrap(~statistic,ncol=4)
                    p <- p+labs(title=paste0("refuges: ",rf,", shortscale: ",ss,", Glac front: ",gl,", Marginal decrease: ",mrg))
                    print(p)
                }
            }
dev.off()

pdf("plot-by-statistic.pdf")
for (stat in unique(dat$statistic))
    for (gl in frnt)
        for (mrg in mrgs)
                {
                    tdf <- filter(dat,statistic==stat,glacfront==gl,marginal.decrease==mrg)
                    tdf$longmean <- as.factor(tdf$longmean)
                    tdf$longmean <- as.factor(tdf$refuges)
                    if (dim(tdf)[1]>0)
                    {
                        p <- tdf %>% ggplot(aes(x=mix,y=value))
                        p <- p+geom_boxplot(aes(group=mix))
                        p <- p + stat_summary(fun.y=median,geom="line",aes(group=1))
                        p <- p + facet_wrap(~longmean+shortscale,ncol=4,labeller="label_both")
                        p <- p+labs(title=paste0("Statistic: ",stat))
                        print(p)
                    }
                }
dev.off()
