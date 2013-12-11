probAnc<-function(p,q,R,mode="budding",analysis="directDesc",Mmax=85,nrep=10000){
	#see Foote, 1996	
	#calculates prob of taxa with indirect desc under budding speciation
		#under infinite time, with p=q or p<q
	#When mode=anagenesis, p is taken to be the rate of anagenesis
	#unused equations related to estimating indirect ancestry? DWB 12-09-13
		#Pinf<-function(p,q,Pp){
		#	res<-numeric()
		#	for(M in 1:80){
		#	res[M]<-Qm(p,q,M)*(1-(1-Pp)^M)
		#		}
		#	sum(res)
		#	}
	#test
	if(!any(mode==c("budding","bifurcating","anagenesis"))){
		stop("Mode not designated, must be 'budding', 'bifurcating' or 'anagenesis'")}
	if(mode=="anagenesis"){message("p will be treated as the rate of anagenesis/pseudospeciation")}
	if(!any(analysis==c("directDesc","indirectDesc"))){
		stop("Analysis type not designated, must be 'directDesc' or 'indirectDesc'")}
	if(nrep<0){stop("Error: Nrep is less than zero?")}
	if(analysis=="directDesc"){
		#get completeness
		Pp<-qsProb2Comp(R=R,p=p,q=q,mode=mode)
		#functions dependent on mode
		if(mode=="budding"){
			Pd<-function(p,q,Ti){exp(-q*(Ti-1))-exp(-q*Ti)}
			PN<-function(p,q,Ti,Ni){(exp(-p*Ti)*((p*Ti)^Ni))/factorial(Ni)}
			#should approximate 2*(p/q)*Pp when R and p are both <<1
			approx<-2*(p/q)*Pp
			maxN<-100
			}
		if(mode=="bifurcating"){
			Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
			PN<-function(p,q,Ti,Ni){
				if(Ni==0){res<-q/(p+q)}
				if(Ni==2){res<-p/(p+q)}
				if(Ni!=2 & Ni!=0){res<-0}
				return(res)
				}
			approx<-2*(p/(p+q))*Pp
			maxN<-2
			}
		if(mode=="anagenesis"){
			Pd<-function(p,q,Ti){exp(-(p+q)*(Ti-1))-exp(-(p+q)*Ti)}
			PN<-function(p,q,Ti,Ni){
				if(Ni==0){res<-q/(p+q)}
				if(Ni==1){res<-p/(p+q)}
				if(Ni!=1 & Ni!=0){res<-0}
				return(res)
				}
			approx<-2*(p/(p+q))*Pp
			maxN<-1
			}
		#now get PA, the probability of a direct descendant of a taxon being sampled
		res<-numeric()
		for(t in 1:nrep){
			Nres<-numeric()
			for(N in 0:maxN){
				Nres[N]<-PN(p=p,q=q,Ti=t,Ni=N)*(1-((1-Pp)^N))
				}
			res[t]<-(1-((1-R)^t))*Pd(p=p,q=q,Ti=t)/Pp*sum(Nres)
			}
		}
	#now prob of sampling an indirect descendant over indefinite time 
		#(nrep in the case of anagenesis)
	if(analysis=="indirectDesc"){
		if(mode=="budding" | mode=="bifurcating"){
			if(p>q){stop(
				"Error: Indirect Descendant formulae are unsolved if p>q, see Foote 1996")}
			Qm<-function(p,q,M){
				x<-(4*p*q)/((p+q)^2)
				res<-((p+q)/(2*p))*(factorial(2*M)/((2^(2*M))*factorial(M)^2))*((x^M)/((2*M)-1))
				return(res)
				}
			Pp<-qsProb2Comp(R=R,q=q,mode="budding")	#as per instructions in Foote, 1996, use 'budding' for both
			#get prob using P'		
			res<-numeric()
			for(M in 1:Mmax){	#85 is the maximum number we can do factorial(2*M) too... but doesn't matter mostly
				res[M]<-Qm(p=p,q=q,M=M)*(1-(1-Pp)^M)
				}
			}
		if(mode=="anagenesis"){
			Pp<-qsProb2Comp(R=R,q=q,mode="anagenesis")
			QmStar<-function(p,q,M,T){
				firstTerm<-numeric()
				for(t in 1:(T-1)){
					firstTerm<-(exp(-q*(t-1))-exp(-q*t))*(exp(-p*t)*((p*t)^(M-1)))/factorial(M-1)
					}
				secondTerm<-exp(-q*(T-1))*(exp(-p*T)*((p*T)^(M-1)))/factorial(M-1)
				res<-sum(firstTerm)+secondTerm
				return(res)
				}
			res<-numeric()
			for(T in 1:nrep){
				Tres<-numeric()
				for(M in 1:Mmax){	
					Tres[M]<-QmStar(p=p,q=q,M=M,T=T)*(1-(1-Pp)^M)
					}
				res[T]<-sum(Tres)
				}
			}
		}
	res<-sum(res)
	names(res)<-NULL
	return(res)
	}
