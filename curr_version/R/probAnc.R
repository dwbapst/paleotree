probAnc<-function(p,q,R){	
	#calculates prob of taxa with indirect desc under budding speciation
		#under infinite time, with p=q or p<q
	Pd<-function(q,Ti){exp(-q*(Ti-1))-exp(-q*Ti)}
	PN<-function(p,Ti,N){(exp(-p*Ti)*(p*Ti)^N)/factorial(N)}
	Qm<-function(p,q,M){
		x<-(4*p*q)/((p+q)^2)
		((p+q)/(2*p))*(factorial(2*M)/((2^(2*M))*factorial(M)^2))*((x^M)/((2*M)-1))
		}
	Pinf<-function(p,q,Pp){
		res<-numeric()
		for(M in 1:80){
		res[M]<-Qm(p,q,M)*(1-(1-Pp)^M)
			}
		sum(res)
		}
	Pp<-qsProb2Comp(R,q)
	res<-numeric()
	for(t in 1:2000){
		Nres<-numeric()
		for(N in 1:100){
			Nres[N]<-PN(p,t,N)*(1-(Pp)^N)
			}
		res[t]<-(((1-(1-R)^t)*Pd(q,t))/Pp)*sum(Nres)
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}