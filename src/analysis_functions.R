
#locus1 is the maternally inherited marker
#all markers have same mutation rate and starting allele freq

#load packages
library("ade4"); library("apex"); library("adegenet"); library("hierfstat");library("pegas"); library("strataG"); library("PopGenReport")

source("make-landscape.R")
source("segment-regression.R")  

#gi
#load("examp-genind.rda")
#tmp
#load("rep.rda")
#crd
#load("popcoord.rda")

#for now HARD CODE..
analysis.func = function(simul_out,n_smp=24,subsample=F,num_loci=10,...){
 
pdf(file="graphics.pdf",width=11,height=5)
par(mfrow=c(2,4))

	#for now, pull out the microsatellite dataset
	gen_ind_obj<-simul_out[[1]]
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
	#calculate all pairwise geographic distances
	pw_geog_dist<-as.matrix(dist(crd[,2:3],upper=T,diag=T))
	#calculate geographic distance FROM ORIGIN (1,1) for our simulations
	dist_origin<-pw_geog_dist[,1]
	#row each population belongs to- row 1 lowest latitude, row 10 highest latitude sort(rep(1:10,10))
	rows_pops<-crd[,2] 
	
	#NUMBER OF ALLELES, across loci, by pop, using adegenet
	sum_stats_gi<-summary(gi_sub)
	#'x' graph
	plot(sum_stats_gi$pop.n.all,rows_pops,ylab="row of landscape",xlab="number of alleles per population")
	#linear graph
	plot(sum_stats_gi$pop.n.all,dist_origin,ylab="distance from origin",xlab="number of alleles per population")
	all_on_dist<-lm(dist_origin~sum_stats_gi$pop.n.all)
	abline(all_on_dist)
	p_val_all<-summary(all_on_dist)[4]$coefficients[8]
	r_sq_all<-as.numeric(summary(all_on_dist)[8])
	alleles_data<-c(as.numeric(coef(all_on_dist)),p_val_all,r_sq_all,var(sum_stats_gi$pop.n.all[1:10]),var(sum_stats_gi$pop.n.all[50:60]),var(sum_stats_gi$pop.n.all[90:100]))
	 names(alleles_data)<-c("intercept","slope","pval","rsq","var_low_lat","var_mid_lat","var_high_lat")

	#HETEROZYG EXPECTED, across loci, by pop using adegenet..
	 temp<-lapply(seppop(gi_sub),summary)
	 het_by_pop<-as.vector(as.numeric(lapply(temp, function(x) mean(x$Hexp))))
	 #lapply(temp, function(x) var(x$Hexp))
	 #'x' graph
	 plot(het_by_pop,rows_pops,ylab="row of landscape",xlab="heterozygosity per population")
	 #linear graph
	 plot(het_by_pop,dist_origin,ylab="distance from origin",xlab="heterozygosity per population")
	 het_on_dist<-lm(dist_origin~het_by_pop)
	 abline(het_on_dist)
	 p_val_het<-summary(het_on_dist)[4]$coefficients[8]
	 r_sq_het<-as.numeric(summary(het_on_dist)[8])
	 het_data<-c(as.numeric(coef(het_on_dist)),p_val_het,r_sq_het,mean(het_by_pop[1:10]),mean(het_by_pop[30:40]),mean(het_by_pop[90:100]))
	 names(het_data)<-c("intercept","slope","pval","rsq","mean_low_lat","mean_mid_lat","mean_high_lat")

	#FAST FST
	all_pw_FST<-pairwise.fstb(gi_sub)
	fst_per_pop_gi<-colMeans(all_pw_FST)
	diag(all_pw_FST)<-NA
	plot(fst_per_pop_gi,rows_pops)
	#color code populations by row!
	#note these are NOT isolation by distance, this is per population FST and distance from origin
	fst_on_dist<-lm(log(fst_per_pop_gi)~log(dist_origin+0.001))
	plot((log(dist_origin)),log(fst_per_pop_gi),pch='')
	pop.num<-1:100
	text((log(dist_origin)),log(fst_per_pop_gi),labels=as.character(pop.num))
	#abline(fst_on_dist)
	#resid(fst_on_dist)
	
	#isolation by distance (IBD)	
	IBD<-lm(c(all_pw_FST)~c(log(pw_geog_dist+0.001)))
	p_val_fst<-summary(IBD)[4]$coefficients[8]
	r_sq_fst<-as.numeric(summary(fst_on_dist)[8])
	plot(pw_geog_dist, all_pw_FST)
	
	#broken stick
	two_reg_stats<-segmentGLM(c(all_pw_FST),c(log(pw_geog_dist+0.001)))
	
	#global FST fitting distribution
	gi_sub_pegas<-genind2loci(gi_sub)
	glob_fst_by_loc<-Fst(gi_sub_pegas)[,2]
	
	#Variance in FST
	var_fst_pop<-as.vector(as.numeric(lapply(FSTpop, function(x) sd(x[,2]))))
	plot(dist_origin,var_fst_pop)
	summary(lm(var_fst_pop~dist_origin))
	abline(lm(var_fst_pop~dist_origin))
	#oddly, variance decreases with distance from origin
	
	fst_data<-c(as.numeric(coef(all_on_dist)),p_val_fst,r_sq_fst,mean(var_fst_pop[1:10]),mean(var_fst_pop[40:50]), mean(var_fst_pop[90:100]))
	names(fst_data)<-c("intercept","slope","pval","rsq","var_low_lat","var_mid_lat","var_high_lat")
	
	dev.off()
	list("ALLELE STATS" = alleles_data, "HETEROZYGOSITY STATS" = het_data, "FST STATS" = fst_data, "TWO REG MODEL" = two_reg_stats, "FIT DISTRIBUTION WEIBULL" = fit.weibull(glob_fst_by_loc), "FIT DISTRIBUTION GAMMA" = fit.gamma(glob_fst_by_loc))
	
}
