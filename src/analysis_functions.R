
#locus1 is the maternally inherited marker
#all markers have same mutation rate and starting allele freq

#load packages
library("ade4"); library("adegenet"); library("hierfstat"); library("pegas"); 	#all needed??
library("strataG"); 	#strataG for mratio and... ??
library("PopGenReport")	#this is the one with fast FST

setwd("C:/Users/shoban/Documents/GitHub/holoSim/src")

source("neighbor_funcs.R")
source("segment-regression.R")  
#source("make-landscape.R")

#tmp is object that contains codominant microsat data in [[1]]
load("../rep.rda")
gi<-tmp[[1]]

#population coordinates crd (pop number, row, column) assumes 100 populations
crd<-matrix(nrow=100,ncol=3)
crd[,1]<-1:100;  crd[,2]<-rep(1:10,each=10);  crd[,3]<-rep(1:10,10)

#ROWS TO FOCUS ON to calculate mean and variance for a statistic, for populations in that row, e.g. bottom row, middle, top row
rows_of_focus<-list(1:10,51:60,91:100)

#calculate all PAIRWISE GEOGRAPHIC DISTANCES among populations- first get rows, columns (coordinates)
pw_geog_dist<-as.matrix(dist(crd[,2:3],upper=T,diag=T))

#calculate geographic DISTANCE FROM ORIGIN (1,1) for our simulations
#this ASSUMES ORIGIN AT 1,1
dist_origin<-pw_geog_dist[,1]

#row each population belongs to- row 1 lowest latitude, row 10 highest latitude 
rows_pops<-crd[,2] 

#names for list of general statistics to calculate
names_stats<-c("intercept","slope","pval","rsq","mean_low_lat","mean_mid_lat","mean_high_lat","var_low_lat","var_mid_lat","var_high_lat")

#for now HARD CODE.. 
analysis.func = function(simul_out,n_smp=24,subsample=F,num_loci=10,doplot=F,...){
 
	if (doplot)   { pdf(file="graphics.pdf",width=11,height=5;  par(mfrow=c(2,4)) }

	#for now, pull out the microsatellite dataset
	gen_ind_obj<-simul_out[[1]]  #probably should use a named object, just in case
	gi<-gen_ind_obj
	
	#NUMBER INDIVIDUALS TO SAMPLE 
	if (!exists("n_smp")) num_samp<-24
	#subset to num_samp individuals
	if (subsample==T) {
		gi_sub<-seppop(gi)
		gi_sub<-lapply(gi_sub, function(x) x[sample(1:nrow(x$tab), num_samp)])
		gi_sub<-repool(gi_sub)
	} else gi_sub<-gi
	gi_sub_gtype<-genind2gtypes(gi_sub)
	
	#TO DO- check which populations have individuals and put in NAS for those with pop size of 0
	
	
	###########  CALCULATING STATISTICS  #############
	
	#NUMBER OF ALLELES, across loci, by pop, using adegenet
		all_by_pop<-summary(gi_sub)$pop.n.all
		alleles_data<-c(lm_summary(dist_origin,all_by_pop), mean_var_rows(all_by_pop))
		names(alleles_data)<-names_stats

	
	#HETEROZYG EXPECTED, across loci, by pop using adegenet..
		temp<-lapply(seppop(gi_sub),summary)
		het_by_pop<-as.vector(as.numeric(lapply(temp, function(x) mean(x$Hexp))))
		het_data<-c(lm_summary(dist_origin,het_by_pop), mean_var_rows(het_by_pop))
		names(het_data)<-names_stats
	 
	 
	#M RATIO, at each locus from strataG
		mrat_by_pop<-colMeans(mRatio(gi_sub_gtype))
		mrat_data<-c(lm_summary(dist_origin,mrat_by_pop), mean_var_rows(het_by_pop))
		names(het_data)<-names_stats
	
	#FAST FST
	#note these are NOT isolation by distance (which is below), this is per population FST and distance from origin
		all_pw_FST<-pairwise.fstb(gi_sub)
		fst_by_pop<-colMeans(all_pw_FST)
		diag(all_pw_FST)<-NA
		fst_data<-c(lm_summary(dist_origin,fst_by_pop), mean_var_rows(fst_by_pop))
		names(fst_by_pop)<-names_stats
	
	#isolation by distance (IBD)	
		diag(pw_geog_dist)<-NA
		IBD_data<-lm_summary(c(pw_geog_dist), c(all_pw_FST))
			
	#broken stick
		two_reg_stats<-segmentGLM(c(pw_geog_dist),log(c(all_pw_FST)))
	
	#global FST fitting distribution
	gi_sub_pegas<-genind2loci(gi_sub)
	glob_fst_by_loc<-Fst(gi_sub_pegas)[,2]
	
	#Variance in FST across loci- does it tell us much?
	all_pw_by_loc<-lapply(seploc(gi_sub), function(x) pairwise.fstb(x))
	concat_pw<-do.call(cbind,lapply(all_pw_by_loc, c))
	var_pw_fst<-apply(concat_pw,1,var)
	#var_fst_pop<-as.vector(as.numeric(lapply(FSTpop, function(x) sd(x[,2]))))
	#plot(dist_origin,var_fst_pop)
	#summary(lm(var_fst_pop~dist_origin)); abline(lm(var_fst_pop~dist_origin))
	#oddly, variance decreases with distance from origin- huh!!!
	#mean(var_fst_pop[1:10]),mean(var_fst_pop[40:50]), mean(var_fst_pop[90:100])
	if (doplot==T) plot(c(pw_geog_dist),var_pw_fst)
#	var_inc<-lm(var_pw_fst~c(pw_geog_dist)) #this takes a long time and is not used later
	#p_val_fst_var<-summary(var_inc)[4]$coefficients[8]
	#r_sq_fst_var<-as.numeric(summary(var_inc)[8])
	
	fst_data<-c(as.numeric(coef(all_on_dist)),p_val_fst,r_sq_fst)
	names(fst_data)<-c("intercept","slope","pval","rsq")
	
	#Nearest Neighbor (FST of the populations nearest to you)
	nn<-near_neighb(all_pw_FST,pw_geog_dist)
	#plot against distance from origin
	plot(nn[,3],nn[,2],ylim=c(0,0.15))
	#plot against population number
	#plot(nn[,1],nn[,2],ylim=c(0,0.15))
	#plot against row number
	#plot(rep(1:10,each=36),nn[,2],ylim=c(0,0.15))
	nn_on_dist<-lm(dist_origin~nn[,2])
	p_val_nn<-summary(nn_on_dist)[4]$coefficients[8]
	r_sq_nn<-as.numeric(summary(nn_on_dist)[8])
	nn_data<-c(as.numeric(coef(nn_on_dist)),p_val_nn,r_sq_nn,mean(nn[1:36,2]),mean(nn[109:144,2]),mean(nn[325:360])) #there should be 360 nearest neighbors??
	names(nn_data)<-c("intercept","slope","pval","rsq","mean_low_lat","mean_mid_lat","mean_high_lat")

		#'x' graph then linear graph
		if (doplot==T) {
			plot(all_by_pop,rows_pops,ylab="row of landscape",xlab="number of alleles per population")
			plot(all_by_pop,dist_origin,ylab="distance from origin",xlab="number of alleles per population");  abline(lm(dist_origin~all_by_pop)) 
			plot(het_by_pop,rows_pops,ylab="row of landscape",xlab="heterozygosity per population")
			plot(het_by_pop,dist_origin,ylab="distance from origin",xlab="heterozygosity per population");  abline(lm(dist_origin~het_by_pop)) 
			plot(fst_by_pop,rows_pops)
			plot(dist_origin,fst_by_pop,pch='');  abline(lm(dist_origin~fst_by_pop))	#plot distance from origin
			pop.num<-1:100;  text(dist_origin,fst_by_pop,labels=as.character(pop.num))
			
			fst_on_dist<-lm_summary(dist_origin,fst_by_pop)	#removed the log here- not needed?
			resid(fst_on_dist)	#anything we can do with the residuals?
		
			plot(pw_geog_dist, all_pw_FST)	#plot IBD
					
		}
		
		
    if (doplot) dev.off()
    list("ALLELE_STATS" = alleles_data,
         "HETEROZYGOSITY_STATS" = het_data,
         "FST_STATS" = fst_data,
         "TWO_REG_MODEL" = two_reg_stats,
         "FIT_DISTRIBUTION_WEIBULL" = fit.weibull(glob_fst_by_loc),
         "FIT_DISTRIBUTION_GAMMA" = fit.gamma(glob_fst_by_loc)
         )
	
}

#Generalized function to extract stats from summaries of the linear models
lm_summary<-function(dist_vector,stat_vector){
	stat_on_dist<-lm(dist_vector~stat_vector)
	#summarize linear model
	coeffic<-(as.numeric(coef(stat_on_dist)))
	p_val<-summary(stat_on_dist)[4]$coefficients[8]
	r_sq<-as.numeric(summary(stat_on_dist)[8])
	c(coeffic, p_val, r_sq)
}

#Generalized function to extract means/ variances from top, middle, bottom rows
mean_var_rows<-function(stats_vector){
	c(
	mean(stats_vector[1:10]),  # var and mean by rows... these three lines need generalizing
    mean(stats_vector[50:60]), #if the landscape is different size
    mean(stats_vector[90:100]), #they will not work
	var(stats_vector[1:10]),  
    var(stats_vector[50:60]), 
    var(stats_vector[90:100])
	)

}