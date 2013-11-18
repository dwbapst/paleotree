pqr2Ps<-function(p,q,r,useExact=TRUE){
	#the probability of sampling at least once an extinct clade of unknown size
	if(useExact){
		#test emily's alternative derivation 10-03-13
		#quadratic solution with minus
		res<-((p+q+r)^2)-(4*p*q)
		1-(((p+q+r)-(sqrt(res)))/(2*p))
	}else{
		#original, needs to be iterated over N
		res<-numeric()
		for(N in 1:1000){
			res1<-((p^(N-1))*(q^N)*choose((2*N)-2,N-1))/(N*((p+q+r)^((2*N)-1)))
			if(is.na(res1)){break}	
			if(res1==Inf){break}
			if(res1==0){break}
			res[N]<-res1
			}
		res<-1-sum(res)
		}
	return(res)
	}