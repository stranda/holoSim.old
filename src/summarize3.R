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

#drop out invariant columns
dat = dat[,-1*which(names(dat) %in% c("shortscale","shortshape","nmut","cpmut","refsize","samp.per.pop", "numloci","computation","dens.scale","nloc","marginal.decrease"))]

treats = unique(dat[,c("refuges","glacfront","longmean","mix")])

mix = do.call(rbind,lapply(1:length(stats),function(i){
    tmpdat = dat[dat$statistic==stats[i],]
    fit = lm(value~mix,data=tmpdat)
    sfit = summary(fit)
    df = t(data.frame(sfit$coefficients[2,c(1,3,4)]))
    rownames(df) = stats[i]
    df
}))


longmean = do.call(rbind,lapply(1:length(stats),function(i){
    tmpdat = dat[dat$statistic==stats[i],]
    fit = lm(value~longmean,data=tmpdat)
    sfit = summary(fit)
    df = t(data.frame(sfit$coefficients[2,c(1,3,4)]))
    rownames(df) = stats[i]
    df
}))

refuges = do.call(rbind,lapply(1:length(stats),function(i){
    tmpdat = dat[dat$statistic==stats[i],]
    fit = aov(value~as.factor(refuges),data=tmpdat)
    sfit = summary(fit)
    df = (data.frame(sfit[[1]][1,c(4,5)]))
    rownames(df) = stats[i]
    df
}))

