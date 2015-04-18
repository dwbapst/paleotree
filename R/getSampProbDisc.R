#' Fit Models of Sampling Probability to Discrete-Interval Taxon Ranges
#' 
#' Uses maximum likelihood to estimate per-interval sampling probability and
#' extinction rate, given a set of discrete-interval taxon ranges from the
#' fossil record. This function can fit models where there are different
#' groupings of taxa with different parameters and different free-floating time
#' intervals with different parameters.
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
#' getSampProbDisc allows for a considerable level of versatility in terms of
#' the variation allowed among taxa in model parameters (extinction rate and
#' sampling probability/rate). Essentially, taxa are divided into different
#' (possibly overlapping) classes which have 'average' parameter values. These
#' average parameters are multiplicatively-combined to calculate per-taxon
#' parameters. For example, perhaps a user hypotheses that taxa which live in a
#' particular environment have different characteristic sampling probabilities,
#' that taxa of several different major clades have different characteristic
#' sampling probabilities and that there may be several temporal shifts in the
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
#' characteristic sampling probability.
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
#' The intervals in timeData should be non-overlapping sequential intervals of
#' roughly equal length. These should be in relative time, so the earliest
#' interval should be 1 and it should increase as the intervals go up with age.
#' This is so differences in interval numbers represents the same rough
#' difference in interval timing. For example, a dataset where all taxa are
#' listed from a set of sequential intervals of similar length, such as North
#' American Mammal assemblage zones, microfossil faunal zones or graptolite
#' biozones can be given as long as they are correctly numbered in sequential
#' order in the input. As a counter example, a dataset which includes taxa
#' resolved only to intervals as wide as the whole Jurassic and taxa resolved
#' to biozones within the Jurassic should not be included in the same input.
#' Drop the taxa from less poorly resolved intervals from such datasets if you
#' want to apply this function, as long as this retains a large enough sample
#' of taxa from the sequential intervals. Note that taxicDivDisc and the "bin_"
#' timescaling methods do not require that intervals be truly sequential (they
#' can be overlapping; see their helpfiles). The output from binTimeData is
#' always sequential, at least by default.
#' 
#' Please check the $message and $convergence elements of the output to make
#' sure that convergence occurred. The likelihood surface can be very flat in
#' some cases, particularly for small datasets (<100 taxa). If convergence was
#' not reached, a warning message is communicated. If the optimizer does not
#' converge, consider increasing iterations or changing the starting values.
#' 
#' @param timeData A 2 column matrix with the first and last occurrences of taxa
#' given in relative time intervals (i.e. ordered from first to last). If a list
#' of length two is given for timeData, such as would be expected if the output 
#' of binTimeData was directly input, the second element is used.

#' @param n_tbins Number of time bins with different sampling/extinction
#' parameters
#' @param grp1 A vector of integers or characters, the same length as the
#' number of taxa in timeData, where each taxon-wise element gives the group ID
#' of the taxon for the respective row of timeData
#' @param grp2 A vector of integers or characters, the same length as the
#' number of taxa in timeData, where each taxon-wise element gives the group ID
#' of the taxon for the respective row of timeData
#' @param est_only If true, function will give back a matrix of ML extinction
#' rates and sampling probabilities per species rather than usual output (see
#' below)
#' @param iterations Maximum number of iterations the optimizer is run for
#' @param initial Values used as initial parameter value to initiate the
#' optimizing function. The same starting value is used for all parameters
#' @return If est_only = TRUE, a matrix of per-taxon sampling and extinction
#' parameters is output.
#' 
#' If est_only = FALSE (default), then the output is a list:
#' 
#' \item{Title}{Gives details of the analysis, such as the number and type of
#' parameters included and the number of taxa analyzed}
#' \item{parameters}{Maximum-likelihood parameters of the sampling model}
#' \item{log.likelihood}{The maximum support (log-likelihood) value}
#' \item{AICc}{The second-order Akaike's Information Criterion, corrected for
#' small sample sizes} \item{convergence}{A number indicating status of
#' convergence; if 0 then convergence was reached; see help file for optim for
#' the respective error list} \item{message}{Messages output by optim; check to
#' make sure that model convergence occurred}
#' 
#' If multi-class models are fit, the element $parameters will not be present,
#' instead there will be several different elements that describe the
#' characteristic parameter 'components' for each class, rather than
#' representing the parameters of actual taxa in that class. As noted in the
#' $Title, these should not be interpreted as the actual rates/probabilities
#' associated with any real taxa but rather as factors which must be multiplied
#' in combination with the estimates for other classes to be meaningful. For
#' example, for taxa of a given group in a given time bin, their extinction
#' rate is the extinction rate component of that time bin times the extinction
#' rate component of their group. Completeness estimates will be output with
#' these parameters as long as classes are not overlapping, as those estimates
#' would not otherwise refer to meaningful groups of taxa.
#' 
#' As the limited time-window option of Foote (1997) is not implemented, any
#' taxa listed as being in a bin with start time 0 and end time 0 (and thus
#' being extant without question) are dropped before the model fitting it
#' performed.

#' @author David W. Bapst, with considerable advice from Michael Foote.

#' @seealso 
#' See the newer version of this method at \code{\link{durationFreq}}.
#'
#' Also see \code{\link{freqRat}}, \code{\link{getSampRateCont}},
#' \code{\link{sProb2sRate}}, \code{\link{qsProb2Comp}}

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
#' #Now let's use binTimeData to bin in intervals of 1 time unit
#' rangesDisc <- binTimeData(rangesCont,int.length=1)
#' 
#' #now, get an estimate of the sampling rate (we set it to 0.5 above)
#' #for discrete data we can estimate the sampling probability per interval (R)
#'     #i.e. this is not the same thing as the instantaneous sampling rate (r)
#' #can use sRate2sProb to see what we would expect
#' sRate2sProb(r=0.5)
#' #expect R = ~0.39
#' 
#' #now we can use maximum likelihood to taxon ranges to get sampling probability
#' SPres1 <- getSampProbDisc(rangesDisc)
#' print(SPres1)      #let's look at the results
#' sProb<-SPres1[[2]][2]
#' print(sProb)
#' #est. R = ~0.42; not too off what we would expect!
#' #for the rate calibrated timescaling methods, we want an estimate of the instantaneous samp rate
#' #we can use sProb2sRate to get the rate. We will also need to also tell it the int.length
#' sRate <- sProb2sRate(sProb,int.length=1)
#' print(sRate)
#' #estimates that r=0.54... Not bad!
#' #Note: for real data, you may need to use an average int.length (no constant length)
#' 
#' \donttest{
#' #this data was simulated under homogenous sampling probabilities, extinction rates
#' #if we fit a model with random groups and allow for multiple timebins
#' 	#AIC should be higher (less informative models)
#' randomgroup <- sample(1:2,nrow(rangesDisc[[2]]),replace=TRUE)
#' SPres2 <- getSampProbDisc(rangesDisc,grp1=randomgroup)
#' SPres3 <- getSampProbDisc(rangesDisc,n_tbins=2)
#' print(c(SPres1$AICc,SPres2$AICc,SPres3$AICc))
#' #and we can see the most simple model has the lowest AICc (most informative model)
#' 
#' #testing temporal change in sampling probabiluty
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=100,maxtaxa=125,maxtime=1000,maxExtant=0,
#'   plot=T)
#' #let's see what the 'true' diversity curve looks like in this case
#' #simulate two sets of ranges at r=0.7 and r=0.1
#' rangesCont <- sampleRanges(taxa,r=1.1)
#' rangesCont2 <- sampleRanges(taxa,r=0.2)
#' #now make it so that taxa which originated after 850 have r=0.1
#' rangesCont[taxa[,3]<850,] <- rangesCont2[taxa[,3]<850,]
#' rangesDisc <- binTimeData(rangesCont)
#' #lets plot the diversity curve
#' taxicDivDisc(rangesDisc)
#' SPres1 <- getSampProbDisc(rangesDisc)
#' SPres2 <- getSampProbDisc(rangesDisc,n_tbins=2)
#' #compare the AICc of the models
#' print(c(SPres1$AICc,SPres2$AICc)) #model 2 looks pretty good
#' #when does it find the break in time intervals?
#' print(rangesDisc[[1]][SPres2$interval.boundaries[2],1])
#' #not so great: estimates 940, not 850 
#'     #but look at the diversity curve: most richness in bin 1 is before 940
#'     #might have found the right break time otherwise...
#'     #the parameter values it found are less great. Finds variation in q	
#' }
#' 
#' @export getSampProbDisc
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
