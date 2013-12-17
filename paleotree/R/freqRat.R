freqRat<-function(timeData,plot=FALSE){
	#timeData is discrete bin data, like from binTimeData
	if(length(timeData)==2){	#if a timeList matrix...
		timeData[[2]][(timeData[[1]][timeData[[2]][,2],2]==0),1]<-NA
		timeData<-timeData[[2]]
		}	
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(timeData,1,diff)<0)){
		stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	durations<-apply(timeData,1,diff)+1
	f1<-sum(durations==1);f2<-sum(durations==2);f3<-sum(durations==3)
	freqRat<-(f2^2)/(f1*f3)
	names(freqRat)<-"freqRat"
	if(plot){hist(durations,breaks=max(durations),xlab="Duration (time-units)",main="")}
	return(freqRat)
	}