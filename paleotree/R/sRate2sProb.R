sRate2sProb<-function(r,int.length=1){
	res<-1-exp(-r*int.length)	#rough estimate
	names(res)<-NULL
	return(res)
	}