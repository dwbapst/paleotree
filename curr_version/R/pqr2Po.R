pqr2Po<-function(p,q,r){
	#the probability of getting an observed event that ends a branch
	res<-1-(2*(1-pqr2Ps(p,q,r))*(p/(p+q)))
	return(res)
	}