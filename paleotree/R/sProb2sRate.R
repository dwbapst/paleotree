sProb2sRate<-function(R,int.length=1){
	res<-(-log(1-R)/int.length)	#rough estimate
	names(res)<-NULL
	return(res)
	}