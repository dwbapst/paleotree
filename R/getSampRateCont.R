#' Fit Models of Sampling Rates to Continuous-Time Taxon Ranges
#' 
#' Uses maximum likelihood to estimate instantaneous sampling and extinction
#' rates, given a set of continuous-time taxon ranges from the fossil record
#' 
#' @details
#' As of version 1.9 of paleotree, this function is retained for historical
#' purposes and users should instead apply the functions listed at
#' \code{\link{durationFreq}}. The new functions do not offer as many options
#' as this function, but are much simpler and (in a sense) offer very different
#' options for constraining parameter space.
#'
#' This function uses maximum-likelihood solutions obtained by Foote (1997).
#' These analyses are ideally applied to data from single stratigraphic section
#' but potentially are applicable to regional or global datasets, although the
#' behavior of those datasets is less well understood.
#' 
#' getSampRateCont allows for a considerable level of versatility in terms of
#' the variation allowed among taxa in model parameters (extinction rate and
#' sampling probability/rate). Essentially, taxa are divided into different
#' (possibly overlapping) classes which have 'average' parameter values. These
#' average parameters are multiplicatively-combined to calculate per-taxon
#' parameters. For example, perhaps a user hypotheses that taxa which live in a
#' particular environment have different characteristic sampling rates, that
#' taxa of several different major clades have different characteristic
#' sampling rates and that there may be several temporal shifts in the
#' characteristic extinction rate. The classification IDs for the first two can
#' be included as per-taxon vectors (of either character or integers) as the
#' arguments 'grp1' and 'grp2' and the hypothesized number of temporal breaks
#' included as the n_tbins argument. A model where taxa differ in parameters
#' across time, clades and environments will thus be fit and the AIC
#' calculated, such that this model's fit can be compared to other (probably
#' simpler) models.
#' 
#' By default, a simple model is fit to the range data where all taxa belong to
#' a single class, with a single characteristic extinction rate and a single
#' characteristic sampling rate.
#' 
#' The timebins option allows for time intervals with different characteristic
#' model parameters which are not defined a priori. The boundaries between time
#' intervals with different characteristic parameters will thus be additional
#' free-floating parameters included in the AIC calculation. If you have the
#' prior inclination that sampling/extinction changed at a specific point in
#' time, then seperate the taxa that originated before and after that point as
#' two different groups and include those classifications as a group in the
#' arguments.
#' 
#' This function does not implement the finite window of observation
#' modification for dealing with range data containing extant taxa and thus
#' that continues up to the recent (Foote, 1997). This is planned for a future
#' update, however. For now, data input into this function should be for taxa
#' that have already gone extinct by the modern and are not presently extant.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' Please check the $message and $convergence elements of the output to make
#' sure that convergence occurred. The likelihood surface can be very flat in
#' some cases, particularly for small datasets (<100 taxa). If convergence was
#' not reached, a warning message is communicated. If the optimizer does not
#' converge, consider increasing iterations or changing the starting values.
#' 
#' As the limited time-window option of Foote (1997) is not implemented, any
#' taxa listed as being at time 0 (and thus being extant) are dropped before
#' the model fitting it performed.
#' 
#' @param timeData Two-column matrix of per-taxon first and last occurrence
#' given in continuous time, relative to the modern (i.e. older dates are also
#' the 'larger' dates).

#' @param n_tbins Number of time bins with different sampling/extinction
#' parameters
#' @param grp1 A vector of integers or characters, the same length as the
#' number of taxa in timeData, where each taxon-wise element gives the group ID
#' of the taxon for the respective row of timeData
#' @param grp2 A vector of integers or characters, the same length as the
#' number of taxa in timeData, where each taxon-wise element gives the group ID
#' of the taxon for the respective row of timeData
#' @param threshold The smallest allowable duration (i.e. the measured difference in
#' the first and last occurrence dates for a given taxon). Durations below this size 
#' will be treated as "one-hit" sampling events.
#' @param est_only If true, function will give back a matrix of ML extinction
#' rates and sampling probabilities per species rather than usual output (see
#' below)
#' @param iterations Maximum number of iterations the optimizer is run for
#' @param initial Values used as initial parameter value to initiate the
#' optimizing function. The same starting value is used for all parameters
#' @return If est_only = T, a matrix of per-taxon sampling and extinction
#' parameters is output.
#' 
#' If est_only = F (default), then the output is a list:
#' 
#' \item{Title}{Gives details of the analysis, such as the number and type of
#' parameters included and the number of taxa analyzed}
#' \item{parameters}{Maximum-likelihood parameters of the sampling model}
#' \item{log.likelihood}{The maximum support (log-likelihood) value}
#' \item{AICc}{The second-order Akaike's Information Criterion, corrected for
#' small sample sizes} \item{convergence}{A number indicating status of
#' convergence; if 0 then convergence was reached; see help file for optim for
#' the respective error list} \item{message}{Messages output by optim(); check
#' to make sure that model convergence occurred}
#' 
#' If multi-class models are fit, the element $parameters will not be present,
#' instead there will be several different elements that describe the
#' characteristic parameter 'components' for each class, rather than
#' representing the parameters of actual taxa in that class. As noted in the
#' $title, these should not be interpreted as the actual rates/probabilities
#' associated with any real taxa but rather as factors which must be multiplied
#' in combination with the estimates for other classes to be meaningful. For
#' example, for taxa of a given group in a given time bin, their extinction
#' rate is the extinction rate component of that time bin times the extinction
#' rate component of their group. Completeness estimates will be output with
#' these parameters as long as classes are not overlapping, as those estimates
#' would not otherwise refer to meaningful groups of taxa.

#' @author David W. Bapst

#' @seealso 
#' See the newer version of this method at \code{\link{durationFreq}}.
#'
#' Also see \code{\link{freqRat}}, \code{\link{getSampProbDisc}},
#' \code{\link{sRate2sProb}}, \code{\link{qsRate2Comp}}

#' @references Foote, M. 1997 Estimating Taxonomic Durations and Preservation
#' Probability. \emph{Paleobiology} \bold{23}(3):278--300.
#' 
#' Foote, M., and D. M. Raup. 1996 Fossil preservation and the stratigraphic
#' ranges of taxa. \emph{Paleobiology} \bold{22}(2):121--140.

#' @examples
#' 
#' #Simulate some fossil ranges with simFossilTaxa
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #now, get an estimate of the sampling rate (we set it to 0.5 above)
#' SRres1 <- getSampRateCont(rangesCont)
#' print(SRres1)   #that's all the results...
#' 
#' #pulling out the sampling rate from the results
#' sRate <- SRres1[[2]][2]
#' print(sRate)	#estimates that sRate=~0.4 (not too bad...)
#' 
#' #this data was simulated under homogenous sampling rates, extinction rates
#' #if we fit a model with random groups and allow for multiple timebins
#' 	#AIC should be higher (less informative)
#' randomgroup <- sample(1:2,nrow(rangesCont),replace=TRUE)
#' SRres2 <- getSampRateCont(rangesCont,grp1=randomgroup)
#' SRres3 <- getSampRateCont(rangesCont,n_tbins=2)
#' SRres4 <- getSampRateCont(rangesCont,n_tbins=3,grp1=randomgroup)
#' print(c(SRres1$AICc,SRres2$AICc,SRres3$AICc,SRres4$AICc))
#' #the most simple model (the first value) has the lowest AICc (most informative!)
#' 
#' @export getSampRateCont
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
