near_neighb<- function(all_pw_FST, pw_geog_dist){
	 nn<-matrix(nrow=8*8*4+8*4*3+4*2,ncol=3)
	 nn_row<-1
	 for (pp in 1:100){
		num_neighb<-length(which((pw_geog_dist[,pp])==1))
		for (n in 1:num_neighb) {
			p1<-as.numeric(which((pw_geog_dist[,pp])==1)[n])
			p2<-pp
			nn[nn_row,2]<-all_pw_FST[p2,p1]; nn[nn_row,1]<-p2
			nn[nn_row,3]<-dist_origin[p2]
			nn_row<-nn_row+1
			}
		}
		nn
}