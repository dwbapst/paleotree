taxicDivCont<-function(timeData,int.length=1,int.times=NULL,plot=TRUE,plotLogRich=FALSE,timelims=NULL,drop.cryptic=FALSE){
	#This function estimates diversity for bins from continuous-time range data
	#input is a per-species matrix of backwards-time FADs and LADs in 2 columns (FADs first)
		#assumes time is in millions of years
	#time interval starts and ends can be pre-input as a 2 column matrix
		#int.length is ignored in this case
	#output (if TRUE) is matrix of bin-start, bit-end, div
	tblen<-int.length
	if(ncol(timeData)==6){	#also allow it to accept taxad objects
		if(!drop.cryptic){
			timeData<-timeData[,3:4,drop=FALSE]
		}else{
			timeDataF<-sapply(unique(timeData[,6]),function(x) max(timeData[x==timeData[,6],3]))
			timeDataL<-sapply(unique(timeData[,6]),function(x) min(timeData[x==timeData[,6],4]))
			timeData<-cbind(timeDataF,timeDataL)
			}
		}	
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	FAD<-as.numeric(timeData[,1]);LAD<-as.numeric(timeData[,2])
	if(is.null(int.times)){
		midtimes<-seq(max(FAD)+2*tblen,min(LAD)-2*tblen,by=-tblen)
		midtimes<-midtimes[midtimes>0]
		int.start<-midtimes+(tblen/2)
		int.end<-midtimes-(tblen/2)
	}else{
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	div<-sapply(1:length(midtimes),function(x) sum(FAD>=int.end[x])-sum(LAD>int.start[x]))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="taxic Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Taxic Richness")
			}
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}