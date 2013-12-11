qsProb2Comp<-function(R,q,p=NULL,mode="budding",nrep=10000){
	#calculate completeness given R and mu
	#based on equations in appendix of Foote, 1996
	if(mode=="budding"){
		Pd<-function(p,q,Ti){exp(-q*(Ti-1))-exp(-q*Ti)}
		}
	if(mode=="bifurcating"){
		Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
		if(is.null(p)){
			p<-q
			message("Origination rate (p) not given, assuming equal to extinction rate")
			}
		}
	if(mode=="anagenesis"){
		Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
		if(is.null(p)){
			p<-q
			message("Rate of pseudo-speciation / anagenesis (p) not given, assuming equal to extinction rate")
			}
		}
	res<-numeric()
	for(t in 1:nrep){
		res[t]<-(1-((1-R)^t))*Pd(p=p,q=q,Ti=t)
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}