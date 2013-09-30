qsRate2Comp<-function(r,q){
	#calculate completeness given r and mu
	res<-r/(r+q)
	names(res)<-NULL
	return(res)
	}