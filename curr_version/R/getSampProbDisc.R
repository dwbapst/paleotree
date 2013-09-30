getSampProbDisc<-function(timeData,n_tbins=1,grp1=NA,grp2=NA,est_only=FALSE,iterations=10000,initial=0.5){
	#this is the multi-parameter maximum likelihood analysis of binned timeData
		#uses a set of binned-interval timeData (just the by-species first and last intervals matrix) 
			#to fit models of different samp probs and ext rates
			#output from binTimeData() can be input directly (only looks at second matrix)
		#can allow for free-moving time windows and different groups
			#these models can then be compared with AIC
		#if est_only=TRUE, then the q and R estimates will be given back per-species
		#this actually runs very slowly; the current settings are optimized for speed (throttle=1). 
			#Increase throttle to 2-4 for inc accuracy
	#x<-runif(100);x<-cbind(x+rexp(100),x);y<-binTimeData(x);getSampProbDisc(y[[2]],est_only=FALSE)
	#x<-runif(100);x<-cbind(x+rexp(100),x);timeData<-binTimeData(x)[[2]];n_tbins=1;grp1=NA;grp2=NA;est_only=FALSE;throttle=1
	#######################
	#REQUIRED FUNCTIONS BELOW
	#######################
	getPp<-function(R,q){
		#calculate completeness given R and q
		res<-numeric()
		for(Ti in 1:10000){
			res[Ti]<-(1-(1-R)^Ti)*(exp(-q*(Ti-1))-exp(-q*Ti))
			}
		sum(res)
		}
	###########################
	f_dur<-function(q,R,dur){
		TMax<-max(dur)
		Ti<-1:TMax	#unique durations
		N<-sapply(Ti,function(t) sum(t==dur))		#number with those durations
		PDT<-exp(-q*(Ti-1))-exp(-q*Ti)
		Rex<-c(1,rep(2,TMax-1))
		ft<- sapply(Ti,function(t) sum(PDT[t:TMax]*((R^Rex[t])*(Ti[t:TMax]-t+1)*((1-R)^(Ti[t:TMax]-t)))))
		ft<-log(ft/sum(ft))	#divide by sum; same as normalizing by Pp (really? COOL!)
		res<-sum((N*ft)[!is.infinite(ft)])
		return(res)
		}
	####################
	make_ft<-function(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s){
		qR_predict<-qR_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
		function(par){
			qR<-qR_predict(par)
			q<-qR[,1]; R<-qR[,2]/10
			data<-cbind(q,R,dur)
			class<-apply(qR,1,function(x) which(apply(unique(qR),1,function(y) all(x==y))))	#get unique classes of q,R
			f<-sapply(unique(class),function(x) f_dur(q[x==class][1],R[x==class][1],dur[x==class]))
			res<-(-sum(f))
			return(res)
			}
		}
	####################
	qR_predict_multpar<-function(FO,LO,n_tbins,grp1,grp2,g1s,g2s){
	#par consists of a vector, with first n_tbins-1 elements represent bin LENGTHS
		#unless n_tbins==1
		#after that, its a matrix with k rows, first col is q, second is R
			#first comes rows associated with each time bin, 
			#than each group in grp1, than grp2
	#THIS FUNCTION SHOULD OUTPUT A Nx2 MATRIX OF q and R VALUES to the support function
		#MWAHAHAH! functions making functions, possible via the magic of lexical scoping!
	#note that little r is used here as a shorthand for big R (samp prob, not samp rate)
	if(n_tbins>1){	#IF THERE ARE TIME BINS
		if(g1s>0){	#IF THERE ARE GROUPS + TIMEBINS
			if(g2s>0){	#IF THERE IS TWO GROUPS +TIMEBINS
				function(par){
					n_tb<-n_tbins-1
					t_raw<-c(0.5,par[1:n_tb])
					t_prop<-t_raw/sum(t_raw)
					t_bl<-t_prop*(max(FO)-min(LO))
					tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
					mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
					mqrt<-matrix(mqr[1:n_tbins,],,2,byrow=TRUE)
					mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]
					mqrg2<-mqr[(n_tbins+g1s+1):(n_tbins+g1s+g2s),]
					tcat<-sapply(FO,function(x) sum(tbin>=x))
					qrt<-sapply(tcat,function(x) mqrt[x,])
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					qrg2<-sapply(grp2,function(x) mqrg2[x,])
					t((qrt+qrg1+qrg2)/3)
				}
			}else{	#IF THERE IS ONE GROUP +TIMEBINS
				function(par){
					n_tb<-n_tbins-1
					t_raw<-c(0.5,par[1:n_tb])
					t_prop<-t_raw/sum(t_raw)
					t_bl<-t_prop*(max(FO)-min(LO))
					tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
					mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
					mqrt<-matrix(mqr[1:n_tbins,],,2,byrow=TRUE)
					mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]
					tcat<-sapply(FO,function(x) sum(tbin>=x))
					qrt<-sapply(tcat,function(x) mqrt[x,])
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					t((qrt+qrg1)/2)
				}
			}
		}else{	#IF THERE ARE NO GROUPS AND TIMEBINS
			function(par){
				n_tb<-n_tbins-1
				t_raw<-c(0.5,par[1:n_tb])
				t_prop<-t_raw/sum(t_raw)
				t_bl<-t_prop*(max(FO)-min(LO))
				tbin<-c(max(FO),max(FO)-cumsum(t_bl))[1:n_tbins]
				mqrt<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
				tcat<-sapply(FO,function(x) sum(tbin>=x))
				tcount<-sapply(sort(unique(tcat)),function(x) sum(tcat==x))
				#if(all(tcount>100)){t(sapply(tcat,function(x) mqrt[x,]))
				#}else{matrix(0.99,length(FO),2,byrow=TRUE)}
				t(sapply(tcat,function(x) mqrt[x,]))
			}
		}
	}else{	#IF THERE ARE NO TIMEBINS
		if(g1s>0){	#IF THERE ARE GROUPS AND NO TIMEBINS
			if(g2s>0){	#IF THERE IS TWO GROUPS AND NO TIMEBINS
				function(par){
					mqr<-matrix(par,,2,byrow=TRUE)
					mqrg1<-mqr[1:g1s,]
					mqrg2<-mqr[(g1s+1):(g1s+g2s),]
					qrg1<-sapply(grp1,function(x) mqrg1[x,])
					qrg2<-sapply(grp2,function(x) mqrg2[x,])
					t((qrg1+qrg2)/2)
				}
			}else{	#IF THERE IS ONE GROUP AND NO TIMEBINS
				function(par){
					mqr<-matrix(par,,2,byrow=TRUE)
					t(sapply(grp1,function(x) mqr[x,]))
				}
			}
		}else{	#IF THERE ARE NO GROUPS AND NO TIMEBINS (2 param model)
			function(par){
				matrix(par,length(FO),2,byrow=TRUE)
				}	
			}
		}			
	}
	#############################
	if(length(timeData)==2){	#if a timeList matrix...
		timeData[[2]][(timeData[[1]][timeData[[2]][,2],1]==0),1]<-NA	#make all extant taxa NAs to LADs
		timeData<-timeData[[2]]
		}	
	#get rid of any NAs
	if(length(grp1)>1){
		if(length(grp1)!=nrow(timeData)){stop("Error: grp1 is not same length as timeData")}
		grp1<-grp1[!is.na(timeData[,1])]

		}
	if(length(grp2)>1){
		if(length(grp2)!=nrow(timeData)){stop("Error: grp2 is not same length as timeData")}		
		grp2<-grp2[!is.na(timeData[,1])]
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(timeData,1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	dur<-apply(timeData,1,diff)+1
	timeData1<-max(timeData)-timeData+1
	FO<-timeData1[,1];LO<-timeData1[,2]
	#NOW CONSTRUCT PROPER SUPPORT FUNCTION
	#NOTE THAT THE FIRST BIN MUST ALWAYS BE THE FIRST INTERVAL
	if(length(grp1)>1){g1s<-length(unique(grp1))}else{g1s<-0}
	if(length(grp2)>1){g2s<-length(unique(grp2))}else{g2s<-0}
	support_ft<-make_ft(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s)
	npar<-(n_tbins-1)+(ifelse(n_tbins>1,2*n_tbins,0)+ifelse(g1s>0,g1s*2,0)
		+ifelse(g2s>0,g2s*2,0))
	if(npar==0){npar<-2}
	par_lim<-c(0.0001,10)
	par_init<-rep(initial,npar)
	par_min<-rep(par_lim[1],npar)
	par_max<-rep(par_lim[2],npar)
	#TIME PARAMS WILL BE MADE PROPORTIONAL TO EACH OTHER TO BE BIN LENS OF TOTAL INTERVAL
	answer<-optim(par_init,support_ft,method="L-BFGS-B",lower=par_min,
		upper=par_max,control=list(maxit=iterations))		#,trace=1
	#answer
	par<-answer$par
	if(answer$convergence!=0){message("Warning: optimizer did not converge; see help page.")}
	if(est_only){
		qR_predict<-qR_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
		qR<-qR_predict(par)
		colnames(qR)<-c("qMax","RMax")
		qR[,2]<-qR[,2]/10
		res<-qR
	}else{
		mes<-answer$message
		conv<-answer$convergence
		SMax<-(-answer$value)
		n<-length(dur)
		aicc<-(2*npar)-(2*SMax)+((2*npar*(npar+1))/(n-npar-1))
		if((n-npar-1)<1){aicc<-"calc failed, npar>(N-1) !"}
		title<-paste("Analysis with",n_tbins,"time bins and",
			sum(c(g1s>0,g2s>0)),"groupings (",g1s,"and",g2s,
			"States),with",npar,"parameters and",n,"taxa")
		title<-c(title,"Note that parameter output is rate components, NOT the typical rate for a member of that group. See help file.")
		if(n_tbins>1){
			n_tb<-n_tbins-1
			t_raw<-c(0.5,par[1:n_tb])
			t_prop<-t_raw/sum(t_raw)
			t_bl<-t_prop*(max(FO)-min(LO))
			t_ends<-c(max(FO),max(FO)-cumsum(t_bl))
			mqR<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
			mqR[,2]<-mqR[,2]/10
			Pp<-sapply(1:nrow(mqR),function(x) getPp(mqR[x,2],mqR[x,1]))
			mqR<-cbind(mqR,Pp)		
			colnames(mqR)<-c("extRate","sampProb","Completeness")
			mqRt<-mqR[1:n_tbins,]
			if(g1s>0){
			mqRg1<-mqR[(n_tbins+1):(n_tbins+g1s),]	
				if(g2s>0){
					mqRg2<-mqR[(n_tbins+1+g1s):(n_tbins+g1s+g2s),]
					res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqRt[,1:2],pars.grp1=mqRg1[,1:2],
						parts.grp2=mqRg2[,1:2],interval.boundaries=t_ends,interval.lengths=t_bl,
						convergence=conv,message=mes)
				}else{
					res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqRt[,1:2],
						pars.grp1=mqRg1[,1:2],interval.boundaries=t_ends,interval.lengths=t_bl,
						convergence=conv,message=mes)
					}
			}else{
				res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.time=mqRt,
					interval.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
				}	
		}else{
			mqR<-matrix(par,,2,byrow=TRUE)
			mqR[,2]<-mqR[,2]/10
			Pp<-sapply(1:nrow(mqR),function(x) getPp(mqR[x,2],mqR[x,1]))
			mqR<-cbind(mqR,Pp)
			colnames(mqR)<-c("extRate","sampProb","Completeness")
			if(g1s>0){
				mqRg1<-mqR[1:g1s,]
				if(g2s>0){
					mqRg2<-mqR[(g1s+1):(g1s+g2s),]
					res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.grp1=mqRg1[,1:2],pars.grp2=mqRg2[,1:2],
						convergence=conv,message=mes)
				}else{
					res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.grp1=mqRg1,convergence=conv,message=mes)
					}
			}else{
				colnames(mqR)<-NULL
				RMax<-mqR[,2]
				qMax<-mqR[,1]
				Pp<-sapply(1:nrow(mqR),function(x) getPp(mqR[x,2],mqR[x,1])) #From Foote, 1997
				res<-list(Title=title[1],parameters=c("extRate"=qMax,"sampProb"=RMax,"Completeness"=Pp),
					log.likelihood=SMax,AICc=aicc,convergence=conv,message=mes)
				}
		}}
	return(res)
	}