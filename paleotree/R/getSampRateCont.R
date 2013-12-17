getSampRateCont<-function(timeData,n_tbins=1,grp1=NA,grp2=NA,threshold=0.1,est_only=FALSE,iterations=1000000,initial=0.5){
	#this is the multi-parameter maximum likelihood analysis of continuous-time fossil ranges
		#uses a set of timeData (FADs and LADs) to fit models of different samp rates and ext rates
		#can allow for free-moving time windows and different groups
			#these models can then be compared with AIC
		#if est_only=TRUE, then the q and r estimates will be given back per-species
	#x<-runif(100);x<-cbind(x+rexp(100),x);getSampRateCont(x)
	#getSampRateCont(x,grp1=(sample(1:2,100,replace=TRUE)))
	#REQUIRED FUNCTIONS BELOW
	######################
	make_ft<-function(dur,FO,LO,n_tbins,grp1,grp2,g1s,g2s){
	qr_predict<-qr_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
	function(par){
		qr<-qr_predict(par)
		q<-qr[,1]; r<-qr[,2]
		ft<-ifelse(dur==0,q/(r+q),q*r*exp(-q*dur)/(r+q))
		-sum(log(ft))
		}
	}
	####################
	qr_predict_multpar<-function(FO,LO,n_tbins,grp1,grp2,g1s,g2s){
	#par consists of a vector, with first n_tbins-1 elements represent bin LENGTHS
		#unless n_tbins==1
		#after that, its a matrix with k rows, first col is q, second is r
			#first comes rows associated with each time bin, 
			#than each group in grp1, than grp2
	#THIS FUNCTION SHOULD OUTPUT A Nx2 MATRIX OF q and r VALUES to the support function
		#MWAHAHAH! functions making functions, possible via the magic of lexical scoping!
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
	#########################
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	#get rid of any NAs
	if(length(grp1)>1){
		if(length(grp1)!=nrow(timeData)){stop("Error: grp1 is not same length as timeData")}
		grp1<-grp1[!is.na(timeData[,1]) & !(timeData[,2]==0)]
		}
	if(length(grp2)>1){
		if(length(grp2)!=nrow(timeData)){stop("Error: grp2 is not same length as timeData")}		
		grp2<-grp2[!is.na(timeData[,1]) & !(timeData[,2]==0)]
		}
	timeData<-timeData[!is.na(timeData[,1]) & !(timeData[,2]==0),]	#drop living taxa
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	FO<-timeData[,1];LO<-timeData[,2]
	dur<-(-LO)-(-FO)
	#THRESHOLD DETERMINES RANGES TOO SMALL TO BE CONSIDERED NOT ONE-TIMERS
	dur[dur<threshold]<-0
	#NOW THAT DUR IS GOOD, CONSTRUCT PROPER SUPPORT FUNCTION
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
		upper=par_max,control=list(maxit=iterations))
	#answer
	if(answer$convergence!=0){message("Warning: optimizer did not converge; see help page.")}
	par<-answer$par
	names(par)<-NULL
	if(est_only){
		qr_predict<-qr_predict_multpar(FO,LO,n_tbins,grp1,grp2,g1s,g2s)
		qr<-qr_predict(par)
		colnames(qr)<-c("qMax","rMax")
		res<-qr
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
		mqr<-matrix(par[(n_tb+1):(length(par))],,2,byrow=TRUE)
		Pp<-mqr[,2]/(mqr[,1]+mqr[,2])
		mqr<-cbind(mqr,Pp)
		colnames(mqr)<-c("extRate","sampRate","Completeness")
		mqrt<-mqr[1:n_tbins,]
		if(g1s>0){
			mqrg1<-mqr[(n_tbins+1):(n_tbins+g1s),]	
			if(g2s>0){
				mqrg2<-mqr[(n_tbins+1+g1s):(n_tbins+g1s+g2s),]
				res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqrt[,1:2],pars.grp1=mqrg1[,1:2],
					pars.grp2=mqrg2[,1:2],int.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
			}else{
				res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.time=mqrt[,1:2],
					pars.grp1=mqrg1[,1:2],int.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
				}
		}else{
			res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.time=mqrt,
				int.boundaries=t_ends,interval.lengths=t_bl,convergence=conv,message=mes)
			}	
	}else{
		mqr<-matrix(par,,2,byrow=TRUE)
		Pp<-mqr[,2]/(mqr[,1]+mqr[,2])
		mqr<-cbind(mqr,Pp)
		colnames(mqr)<-c("extRate","sampRate","Completeness")
		if(g1s>0){
			mqrg1<-mqr[1:g1s,]
			if(g2s>0){
				mqrg2<-mqr[(g1s+1):(g1s+g2s),]
				res<-list(Title=title,log.likelihood=SMax,AICc=aicc,pars.grp1=mqrg1[,1:2],pars.grp2=mqrg2[,1:2],
					convergence=conv,message=mes)
			}else{
				res<-list(Title=title[1],log.likelihood=SMax,AICc=aicc,pars.grp1=mqrg1,message=mes)
				}
		}else{
			colnames(mqr)<-NULL
			rMax<-mqr[,2]
			qMax<-mqr[,1]
			Pp<-rMax/(rMax+qMax) #From Solow and Smith, 1997
			res<-list(Title=title[1],parameters=c("extRate"=qMax,"sampRate"=rMax,"Completeness"=Pp),
				log.likelihood=-SMax,AICc=aicc,convergence=conv,message=mes)
			}
		}}
	return(res)
	}