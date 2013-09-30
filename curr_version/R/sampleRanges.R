sampleRanges<-function(taxad,r,alpha=1,beta=1,rTimeRatio=1,modern.samp.prob=1,min.taxa=2,
	ranges.only=TRUE,minInt=0.01,merge.cryptic=TRUE,randLiveHat=TRUE,alt.method=FALSE,plot=FALSE){
	#sample ranges using a taxad matrix as input 
	#if (ranges.only=TRUE): outputs matrix of FADs/LADs, with NAs for unsampled taxa
	#if (ranges.only=FALSE): outputs per-species list with vectors of dates where that species was sampled
		#ranges and occurance are output on a BACKWORD-moving timescale as expected for paleo data
	#if modern.samp.prob=1, then all still-living taxa (taxa at 0 for LAD) are ALWAYS last observed at zero
		#this approximates the fact that we think the present-day living biota is almost perfectly sampled
			#(well, relative to the modern)
	#if r is a vector, then it is considered to be a vector of per-species sampling rates
		#r is the mean sampling rate for a species
	#rTimeRatio is the ratio of the latest taxon-avg sampling rates by the earliest taxonavg sampling rates
		#i.e. the proportional increase over (a) the entire taxon's history or (b) the taxon duration (if specific)
		#either a single value for all taxa or taxon-specific values
	#all parameters can be given as single values or species-specific values
	#names<-paste("t",1:4,sep="");taxad<-cbind(c(250,230,210,200),c(240,215,205,0))
	#min.taxa=2;minInt=0.01;modern.samp.prob=0;plot=T;ranges.only=F;alt.method=F;randLiveHat=TRUE;merge.cryptic=TRUE
	#r<-c(0.2,0.1,0.3,0.4);alpha<-4;beta<-4;rTimeRatio<-2
	#r<-c(0,0.1,0.3,0.4);alpha<-beta<-rTimeRatio<-2
	if(ncol(taxad)==6){				#also allow it to accept taxad objects
		living<-taxad[,5]
		cryptic<-sapply(1:nrow(taxad),function(x) if(taxad[x,1]!=taxad[x,6]){which(taxad[,1]==taxad[x,6])}else{NA})
		timeData<-taxad[,3:4,drop=FALSE]
		names<-if(is.null(rownames(taxad))){paste("t",taxad[,1],sep="")}else{rownames(taxad)}
		}
	if(ncol(taxad)==2){			#assumes it has two matrices
		living<-rep(0,nrow(taxad))
		living[taxad[,2]==0]<-1
		cryptic<-rep(NA,nrow(taxad))
		timeData<-taxad
		names<-if(is.null(rownames(taxad))){paste("t",1:nrow(taxad),sep="")}else{rownames(taxad)}
		}
	if(nrow(taxad)<min.taxa){stop("Error: min.taxa set higher than number of taxa in input")}	
	if(merge.cryptic==TRUE & nrow(taxad)<sum(!is.na(cryptic))){
		stop("Error: min.taxa set higher than number of non-cryptic taxa in input")}
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	#check input parameters for inconsistencies
	if(any(r<0)){stop("Error: some r values < 0")}
	if(any(alpha<=0) | any(beta<=0)){stop("Error: some shape parameters (alpha, beta) <= 0")}
	if(length(alpha)==1){alpha<-rep(alpha,length(names))
		}else{if(length(alpha)!=length(names)){
		stop("Error: Multiple alpha values input but not same length as number of taxa?")}}
	if(length(beta)==1){beta<-rep(beta,length(names))
		}else{if(length(beta)!=length(names)){
		stop("Error: Multiple beta values input but not same length as number of taxa?")}}
	if(length(rTimeRatio)!=1){if(length(rTimeRatio)!=length(names)){
		stop("Error: Multiple rTimeRatio values input but not same length as number of taxa?")}}
	midpoint<-(max(timeData)+min(timeData))/2
	taxa.midpoints<-((timeData[,1]+timeData[,2])/2)
	#need to calculate rTimeChange based on rTimeRatio
	r_start<-r/((rTimeRatio/2)+0.5)
	r_end<-2*r*(1-1/(rTimeRatio+1))
	if(length(r)==1 & length(rTimeRatio)==1){
		#thanks to Emily King and Brian Koch for helping me with this bit!
		rTimeChange<-(r_end-r_start)/(max(timeData)-min(timeData))
		rTimeChange<-rep(rTimeChange,length(names))
	}else{	#if there's multiple r values, multiple ratios or both...
		dur1<-timeData[,1]-timeData[,2]
		rTimeChange<-(r_end-r_start)/dur1
		}
	if(length(r)==1){		#get per-taxon r.avg for each taxon, if not input already
		taxa.timeChange<-midpoint-taxa.midpoints
		r<-r+(rTimeChange*taxa.timeChange)
		}
	if(length(r)!=length(names)){stop("Error: Multiple r values given but not same length as number of taxa?")}
	if(any(r<0)){stop("Error: rTimeRatio so low as to cause some species-specific avg r < 0")}
	if(any(alpha!=1) | any(beta!=1) | any(rTimeChange!=0) | alt.method){		#get per species time vectors of r
		rangesTimes<-apply(timeData,1,function(x) seq(x[1],x[2],by=-minInt))
		#get the time change component per minInt of the taxon ranges
		#first get the different of each min int from the midpoints, then calculate the change r accordint to rTimeChange
		taxa.shiftTimes<-lapply(1:length(r),function(x) taxa.midpoints[x]-rangesTimes[[x]])
		taxa.rShiftTime<-lapply(1:length(r),function(x) rTimeChange[x]*taxa.shiftTimes[[x]])
		#if randLiveHat is TRUE, scale 
		#for each taxon decide on an end.hat; use runif if extant, 1 is extinct
		if(randLiveHat){end.hat<-ifelse(living==1,runif(length(living)),1)
			}else{end.hat<-rep(1,length(rangesTimes))}
		#now get the hat component, combine with the time component
		scaledTimes<-lapply(1:length(end.hat),function(x) seq(0,end.hat[x],length.out=length(rangesTimes[[x]])))
		rHat<-lapply(1:length(r),function(x) r[x]*dbeta(scaledTimes[[x]],alpha[x],beta[x]))
		rHatTime<-lapply(1:length(r),function(x) rHat[[x]]+taxa.rShiftTime[[x]])
		#transform so that the hats always are above zero but have correct means and correct timechange slope (so complicated!)
		#rHatTime<-lapply(rHatTime,function(x) (x-min(x))/(1-(min(x)/mean(x))))
		rHatTime1<-lapply(rHatTime,function(x) ifelse(x<0,0,x))
		rHatTime<-lapply(1:length(rHatTime),function(x) if(mean(rHatTime1[[x]])>0){
				mean(rHatTime[[x]])*rHatTime1[[x]]/mean(rHatTime1[[x]])
			}else{mean(rHatTime[[x]])*rHatTime1[[x]]}
			)
	}else{	
		rHatTime<-NULL
		}
	if(plot){if(!is.null(rHatTime)){
		plot(0:1,c(min(unlist(rHatTime)),max(unlist(rHatTime))),type="n",xlim=c(max(timeData),min(timeData)),	
			xlab="Time (Time-Units before Present)",ylab=c("Instant. Sampling Rate","(per lineage time-units)"))
		for(i in 1:length(r)){lines(rangesTimes[[i]],rHatTime[[i]],lwd=2)}
	}else{
		plot(c(max(timeData),min(timeData)),c(min(r),max(r)),type="n",xlim=c(max(timeData),min(timeData)),
			xlab="Time (Time-Units before Present)",ylab=c("Instant. Sampling Rate","(per lineage time-units)"))
		for(i in 1:length(r)){lines(timeData[i,],c(r[i],r[i]),lwd=2)}
		}}
	#plot(rangesTimes[[2]],rHatTime[[2]],type="l",lwd=3,xlim=c(max(rangesTimes[[2]]),
	#	min(rangesTimes[[2]])),xlab="time",ylab="Instantaneous Sampling Rate (per Lmy)")
	#(x<-cbind(r,sapply(rHatTime,mean)));plot(x);abline(0,1)	#means should be close to actual r
	#okay, now to actually do the sampling
	redo<-TRUE
	while(redo){
		samp_occ<-list()	#sampled occurances
		for(i in 1:nrow(timeData)){
			if(is.null(rHatTime)){
				if(r[i]>0){
					samps<-timeData[i,1]	#the time of speciation is the lower bound
					while(min(samps)>timeData[i,2]){	#keep sampling until you go past the time of extinction
						samps<-c(samps,min(samps)-rexp(1,rate=r[i]))}
					samps<-samps[-c(1,length(samps))]
				}else{samps<-numeric(0)}
			}else{
				samps<-rangesTimes[[i]][sapply(rHatTime[[i]],function(x) runif(1)<=(x*minInt))]
				}
			if(modern.samp.prob>0){if(living[i]==1){		#rework modern.sampling as a probability
				if(runif(1)<=modern.samp.prob){samps<-c(samps,0)}
				}}
			if(length(samps)>0){samp_occ[[i]]<-samps}else{samp_occ[[i]]<-NA}
			}
		redo<-min.taxa>sum(sapply(samp_occ,function(x) !any(is.na(x))))
		}
	names(samp_occ)<-names
	if(any(!is.na(cryptic)) & merge.cryptic){for(i in which(!is.na(cryptic))){
		samp_occ[[cryptic[i]]]<-c(samp_occ[[cryptic[i]]],samp_occ[[i]])
		samp_occ[[i]]<-NA
		}}
	if(ranges.only){
		ranges<-cbind(sapply(samp_occ,max),sapply(samp_occ,min))
		rownames(ranges)<-names;colnames(ranges)<-c("FAD","LAD")
		res<-ranges
	}else{
		res<-samp_occ
		}
	return(res)
	}