
#locus1 is the maternally inherited marker
#all markers have same mutation rate and starting allele freq

#load packages
library("ade4"); library("apex"); library("adegenet"); library("hierfstat");library("pegas"); library("strataG"); library("PopGenReport")

#source("make-landscape.R")
#source("segment-regression.R")  

#gi
#load("examp-genind.rda")
#tmp
#load("rep.rda")
#crd
#load("popcoord.rda")

#for now HARD CODE..
analysis.func = function(simul_out,n_smp=24,subsample=F,num_loci=10,doplot=F,...){
 
if (doplot)
    {
        pdf(file="graphics.pdf",width=11,height=5)
        par(mfrow=c(2,4))
    }
	#for now, pull out the microsatellite dataset
	gen_ind_obj<-simul_out[[1]]  #probably should use a named object, just in case
	gi<-gen_ind_obj
	
	#number individuals to sample 
	if (!exists("n_smp")) num_samp<-24
	#subset to num_samp individuals
	if (subsample==T) {
		gi_sub<-seppop(gi)
		gi_sub<-lapply(gi_sub, function(x) x[sample(1:nrow(x$tab), num_samp)])
		gi_sub<-repool(gi_sub)
	} else gi_sub<-gi
	
	#calculate distances among populations- first get rows, columns (coordinates)
	crd<-landscape.popcoord(simul_out[[4]])
	#TO DO- check which populations have individuals and put in NAS for those with pop size of 0
	#calculate all pairwise geographic distances
	pw_geog_dist<-as.matrix(dist(crd[,2:3],upper=T,diag=T))
	#calculate geographic distance FROM ORIGIN (1,1) for our simulations
        #this will need to be generalized for different refugia locations
	dist_origin<-pw_geog_dist[,1]
	#row each population belongs to- row 1 lowest latitude, row 10 highest latitude sort(rep(1:10,10))
	rows_pops<-crd[,2] 
	
	#NUMBER OF ALLELES, across loci, by pop, using adegenet
	sum_stats_gi<-summary(gi_sub)
	#'x' graph
	if (doplot==T) plot(sum_stats_gi$pop.n.all,rows_pops,ylab="row of landscape",xlab="number of alleles per population")
	#linear graph
	if (doplot==T) plot(sum_stats_gi$pop.n.all,dist_origin,ylab="distance from origin",xlab="number of alleles per population")
	all_on_dist<-lm(dist_origin~sum_stats_gi$pop.n.all)
	if (doplot==T) abline(all_on_dist)
	p_val_all<-summary(all_on_dist)[4]$coefficients[8]
	r_sq_all<-as.numeric(summary(all_on_dist)[8])
	alleles_data<-c(as.numeric(coef(all_on_dist)),
                        p_val_all,
                        r_sq_all,
                        var(sum_stats_gi$pop.n.all[1:10]),  #these three lines need generalizing
                        var(sum_stats_gi$pop.n.all[50:60]), #if the landscape is different size
                        var(sum_stats_gi$pop.n.all[90:100]))#they will not work
	 names(alleles_data)<-c("intercept","slope","pval","rsq","var_low_lat","var_mid_lat","var_high_lat")

	#HETEROZYG EXPECTED, across loci, by pop using adegenet..
	 temp<-lapply(seppop(gi_sub),summary)
	 het_by_pop<-as.vector(as.numeric(lapply(temp, function(x) mean(x$Hexp))))
	 #lapply(temp, function(x) var(x$Hexp))
	 #'x' graph
	 if (doplot==T) plot(het_by_pop,rows_pops,ylab="row of landscape",xlab="heterozygosity per population")
	 #linear graph
	 if (doplot==T) plot(het_by_pop,dist_origin,ylab="distance from origin",xlab="heterozygosity per population")
	 het_on_dist<-lm(dist_origin~het_by_pop)
	 if (doplot==T) abline(het_on_dist)
	 p_val_het<-summary(het_on_dist)[4]$coefficients[8]
	 r_sq_het<-as.numeric(summary(het_on_dist)[8])
	 het_data<-c(as.numeric(coef(het_on_dist)),p_val_het,r_sq_het,mean(het_by_pop[1:10]),mean(het_by_pop[30:40]),mean(het_by_pop[90:100]))
	 names(het_data)<-c("intercept","slope","pval","rsq","mean_low_lat","mean_mid_lat","mean_high_lat")

	#FAST FST
	all_pw_FST<-pairwise.fstb(gi_sub)
	fst_per_pop_gi<-colMeans(all_pw_FST)
	diag(all_pw_FST)<-NA
	if (doplot==T) plot(fst_per_pop_gi,rows_pops)
	#color code populations by row!
	#note these are NOT isolation by distance, this is per population FST and distance from origin
	fst_on_dist<-lm(log(fst_per_pop_gi)~log(dist_origin+0.001))
	if (doplot==T) plot((log(dist_origin)),log(fst_per_pop_gi),pch='')
	pop.num<-1:100
	if (doplot==T) text((log(dist_origin)),log(fst_per_pop_gi),labels=as.character(pop.num))
	#abline(fst_on_dist)
	#resid(fst_on_dist)
	
	#isolation by distance (IBD)	
	diag(pw_geog_dist)<-NA
	IBD<-lm(c(all_pw_FST)~c(log(pw_geog_dist+0.001)))
	p_val_fst<-summary(IBD)[4]$coefficients[8]
	r_sq_fst<-as.numeric(summary(IBD)[8])
	if (doplot==T) plot(pw_geog_dist, all_pw_FST)
	
	#broken stick
	two_reg_stats<-segmentGLM(c(pw_geog_dist),log(c(all_pw_FST)))
	
	#global FST fitting distribution
	gi_sub_pegas<-genind2loci(gi_sub)
	glob_fst_by_loc<-Fst(gi_sub_pegas)[,2]
	
	#Variance in FST
	all_pw_by_loc<-lapply(seploc(gi_sub), function(x) pairwise.fstb(x))
	concat_pw<-do.call(cbind,lapply(all_pw_by_loc, c))
	var_pw_fst<-apply(concat_pw,1,var)
	#var_fst_pop<-as.vector(as.numeric(lapply(FSTpop, function(x) sd(x[,2]))))
	#plot(dist_origin,var_fst_pop)
	#summary(lm(var_fst_pop~dist_origin)); abline(lm(var_fst_pop~dist_origin))
	#oddly, variance decreases with distance from origin- oops!
	#mean(var_fst_pop[1:10]),mean(var_fst_pop[40:50]), mean(var_fst_pop[90:100])
	if (doplot==T) plot(c(pw_geog_dist),var_pw_fst)
#	var_inc<-lm(var_pw_fst~c(pw_geog_dist)) #this takes a long time and is not used later
	#p_val_fst_var<-summary(var_inc)[4]$coefficients[8]
	#r_sq_fst_var<-as.numeric(summary(var_inc)[8])
	
	fst_data<-c(as.numeric(coef(all_on_dist)),p_val_fst,r_sq_fst)
	names(fst_data)<-c("intercept","slope","pval","rsq")
	
    if (doplot) dev.off()
    list("ALLELE_STATS" = alleles_data,
         "HETEROZYGOSITY_STATS" = het_data,
         "FST_STATS" = fst_data,
         "TWO_REG_MODEL" = two_reg_stats,
         "FIT_DISTRIBUTION_WEIBULL" = fit.weibull(glob_fst_by_loc),
         "FIT_DISTRIBUTION_GAMMA" = fit.gamma(glob_fst_by_loc)
         )
	
}
