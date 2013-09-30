binTimeData<-function(timeData,int.length=1,start=NA,int.times=NULL){
	#bin temporal data
	#input: continuous time data (two column of FADs and LADs)
	#output: a list with two 2-col matrices as elements, bin-times and taxon occurences
			#intervals, UNLIKE TIME, always go up (earliest is 1 and increase...)
		#arbitrarily starts bin at the first fad; this can be changed by setting 'start'
			#start must be greater than max(timeData)
			#the last bin is cut off at zero (present day)
	#x<-c(0,runif(99));timeData<-cbind(x+rexp(100),x);int.length=1;start=NA;int.times=NULL
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data?")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	if(is.null(int.times)){
		if(is.na(start)){start<-max(timeData)+int.length}else{if(start<max(timeData)){stop("Error:Start<max(timeData)?")}}
		end<-start-(ceiling((start-min(timeData))/int.length)+1)*int.length
		bins<-seq(start,end,by=-int.length)
		bins<-unique(ifelse(bins<0,0,bins))	#get rid of any extra zeroes or negative numbers
		fads<-sapply(timeData[,1],function(x) which(bins<x)[1]-1)
		lads<-sapply(timeData[,2],function(x) which(bins<x)[1]-1)
		if(any(timeData[,1]==0) | any(timeData[,2]==0)){
			bins<-c(bins,0)
			fads[timeData[,1]==0]<-length(bins)-1
			lads[timeData[,2]==0]<-length(bins)-1
			}
		res<-list(int.times=cbind(int.start=bins[1:(length(bins)-1)],int.end=bins[2:length(bins)]),
			taxon.times=cbind(first.int=fads,last.int=lads))
	}else{
		int.durs<-int.times[,1]-int.times[,2]
		if(any(int.durs<=0)){stop("Error: Some input time intervals have zero or negative durations?")}
		int.times<-int.times[order(int.durs),]
		Fint<-sapply(timeData[,1],function(x) which(apply(int.times,1,function(y) y[1]>=x & y[2]<x))[1])
		Lint<-sapply(timeData[,2],function(x) which(apply(int.times,1,function(y) y[1]>=x & y[2]<x))[1])
		if(any(int.times[,2]==0)){
			Fint[timeData[,1]==0]<-which(int.times[,2]==0)[1]
			Lint[timeData[,2]==0]<-which(int.times[,2]==0)[1]
			}
		taxon.times<-cbind(first.int=Fint,last.int=Lint)
		rownames(taxon.times)<-rownames(timeData)
		taxon.times<-taxon.times[!apply(taxon.times,1,function(x) any(is.na(x))),]
		new.order<-rank(-int.times[,1])	
		taxon.times[,1]<-new.order[taxon.times[,1]]
		taxon.times[,2]<-new.order[taxon.times[,2]]
		int.times<-int.times[order(-int.times[,1]),]
		res<-list(int.times=int.times,taxon.times=taxon.times)
		}
	return(res)
	}