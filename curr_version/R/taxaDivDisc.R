taxicDivDisc<-function(timeList,int.times=NULL,plot=TRUE,plotLogRich=FALSE,timelims=NULL,extant.adjust=0.01,split.int=TRUE){
	#this function estimates diversity for binned intervals from discrete interval range data
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
	#time interval starts and ends can be pre-input as a 2 column matrix
		#HOWEVER this could be pretty misleading!
		#standing richness may never be high as the apparent richness of some bins
	#output (if TRUE) is 3 col matrix of int-start, int-end, div
	if(class(timeList[[1]])!="matrix"){if(class(timeList[[1]])=="data.frame"){timeList[[1]]<-as.matrix(timeList[[1]])
		}else{stop("Error: timeList[[1]] not of matrix or data.frame format")}}
	if(class(timeList[[2]])!="matrix"){if(class(timeList[[2]])=="data.frame"){timeList[[2]]<-as.matrix(timeList[[2]])
		}else{stop("Error: timeList[[2]] not of matrix or data.frame format")}}
	intMat<-timeList[[1]]	#the intervals the DATA is given in
	timeData<-timeList[[2]]
	#if(drop.extant){timeData[[2]][(timeData[[1]][timeData[[2]][,2],1]==0),1]<-NA}
	intMat[intMat[,1]==0,1]<-extant.adjust
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(intMat,1,diff)>0)){stop("Error: timeList[[1]] not in intervals in time relative to modern")}
	if(any(intMat[,2]<0)){stop("Error: Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeData,1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	Fint<-as.numeric(timeData[,1]);Lint<-as.numeric(timeData[,2])
	FAD<-intMat[Fint,1];LAD<-intMat[Lint,2]
	if(is.null(int.times)){
		avg_dur<-abs(mean(apply(timeList[[1]],1,diff)))
		int.bounds<-unique(c(intMat,max(intMat)+avg_dur,min(intMat)-avg_dur))		#add a little space at start and end
		int.bounds<-int.bounds[order(-int.bounds)]
		intMat<-cbind(int.bounds[-length(int.bounds)],int.bounds[-1])
		int.start<-intMat[,1];int.end<-intMat[,2]
		midtimes<-apply(intMat,1,mean)
	}else{
		if(split.int){	#if split.int, then any interval times given are split at discrete time intervals
			splinters<-sort(unique(c(intMat)))
			mustSplit<-apply(int.times,1,function(x) any(sapply(splinters,function(y) x[1]>y & x[2]<y)))
			if(any(mustSplit)){
				for(i in which(mustSplit)){
					splitter<-splinters[sapply(splinters,function(y) int.times[i,1]>y & int.times[i,2]<y)]
					for(j in splitter){		#in case there is more than one splitter
						int.times<-rbind(int.times,c(int.times[i,1],j),c(j,int.times[i,2]))
						}
					}
				int.times<-int.times[-which(mustSplit),]
				int.times<-int.times[order(-int.times[,1]),]
				}
			}
		int.start<-int.times[,1];int.end<-int.times[,2]
		midtimes<-(int.start+int.end)/2
		}
	div<-sapply(1:length(midtimes),function(x) sum(FAD>int.end[x])-sum(LAD>=int.start[x]))
	#div<-sapply(min(timeData):max(timeData),function(x) 	sum(FAD<=x & LAD>=x))
	if(plot){
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		div1<-c(div,div)[order(times1)]
		times1<-sort(times1)
		if(plotLogRich){
			plot(times1[div1>0],div1[div1>0],type="l",log="y",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Taxic Richness (Log Scale)")		
		}else{
			plot(times1,div1,type="l",
				xlim=if(is.null(timelims)){c(max(times1),max(0,min(times1)))}else{timelims},
				xlab="Time (Before Present)",ylab="Taxic Richness")
			}
		}
	res<-cbind(int.start,int.end,int.div=div)
	return(invisible(res))
	}