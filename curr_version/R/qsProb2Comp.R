qsProb2Comp<-function(R,q){
	#calculate completeness given R and mu
	res<-numeric()
	for(t in 1:10000){
		res[t]<-(1-(1-R)^t)*(exp(-q*(t-1))-exp(-q*t))
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}