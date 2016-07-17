#These functions are used by the analysis functions to summarize the genetic summary statistics

#Generalized function to extract stats from summaries of the linear models
lm_summary<-function(dist_vector,stat_vector){
	stat_on_dist<-lm(dist_vector~stat_vector)
	#summarize linear model
	coeffic<-(as.numeric(coef(stat_on_dist)))
	p_val<-summary(stat_on_dist)[4]$coefficients[8]
	r_sq<-as.numeric(summary(stat_on_dist)[8])
	c(coeffic, p_val, r_sq)
}

#Generalized function to extract MEANS/ VARIANCES from top, middle, bottom ROWS
#this take rows of focus and applies mean, var to them
#should be ok for different size landscapes and as many rows as you want to look at
mean_var_rows<-function(stats_vector, rows_focus){
	c(
	sapply(rows_of_focus, function (x) mean(all_by_pop[x])),
	sapply(rows_of_focus, function (x) var(all_by_pop[x]))
	)
}

#General function to get summary stats glued together and name the vector
get_sum_stats<-function(stat_vector, dist_vector, rows_list){
	stats_data<-c(lm_summary(dist_vector,stat_vector), mean_var_rows(stat_vector, rows_list))
	names(stats_data)<-names_stats
	stats_data
}
