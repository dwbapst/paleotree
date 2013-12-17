perCapitaRates<-function(timeList,plot=TRUE,plotLogRates=FALSE,drop.extant=FALSE,isExtant=NULL,showLegend=TRUE){
	#
	#this function estimates per-capita rates for binned intervals from discrete interval range data
		#based on Foote, 2000
	#input is a list with (1) interval times matrix and (2) species FOs and LOs
	#time interval starts and ends can be pre-input as a 2 column matrix
		#HOWEVER this could be pretty misleading!
		#standing richness may never be high as the apparent richness of some bins
	#if there is a single interval that is start=0 and end=0, it is dropped and all appearances in this interval are either
		#(a) melded onto the next most recent interval
		#(b) dropped if drop.extant=TRUE
	#output is a matrix of int-start, int-end, p, q
	#example
		#timeList<-rangesDisc;int.times=NULL;plot=TRUE;plotLogRates=FALSE;timelims=NULL
			#drop.extant=FALSE;isExtant=NULL;split.int=TRUE
	#
	if(class(timeList[[1]])!="matrix"){if(class(timeList[[1]])=="data.frame"){timeList[[1]]<-as.matrix(timeList[[1]])
		}else{stop("Error: timeList[[1]] not of matrix or data.frame format")}}
	if(class(timeList[[2]])!="matrix"){if(class(timeList[[2]])=="data.frame"){timeList[[2]]<-as.matrix(timeList[[2]])
		}else{stop("Error: timeList[[2]] not of matrix or data.frame format")}}
	intMat<-timeList[[1]]	#the intervals the DATA is given in
	timeData<-timeList[[2]]
	timeData<-timeData[!is.na(timeData[,1]),,drop=FALSE]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(intMat,1,diff)>0)){stop("Error: timeList[[1]] not in intervals in time relative to modern")}
	if(any(intMat[,2]<0)){stop("Error: Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeData,1,diff)<0)){stop("Error: timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeList[[2]] <0 ?")}
	#get rid of modern intervals
	modInt<-which(intMat[,1]==0)
	if(length(modInt)>0){
		if(drop.extant){
			timeData<-timeData[!(timeData[,2]==modInt),]
		}else{
			isExtant<-timeData[,2]==modInt
			altMod<-which(!intMat[,1]==0 & intMat[,2]==0)
			timeData[timeData[,1]==modInt,1]<-altMod
			timeData[timeData[,2]==modInt,2]<-altMod
			}
		intMat<-intMat[-modInt,]
		}
	if(is.null(isExtant)){isExtant<-rep(FALSE,nrow(timeData))}	
	# make sure data is absolutely sequential and overlapping
	if(!(all(apply(intMat,2,function(x) identical(x,rev(sort(x))) & identical(x,unique(x)))) &
	   all(sapply(2:nrow(intMat),function(x) intMat[x,1]==intMat[x-1,2])))){
		stop("Error: Sorry, intervals need to be sequential, ordered and non-overlapping.")}
	#now get int length to modify rates with
	intlen<-(-apply(intMat,1,diff))
	#modify timeData if there are any extant taxa
	if(any(isExtant)){timeData[isExtant,2]<-nrow(intMat)+1}
	#get the four values
	Nbt<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]<x & timeData[,2]>x))
	NbL<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]<x & timeData[,2]==x))
	NFt<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]==x & timeData[,2]>x))	
	NFL<-sapply(1:nrow(intMat),function(x) sum(timeData[,1]==x & timeData[,2]==x))
	N<-Nbt+NbL+NFt+NFL	#total diversity
	#now calculate rates
	pRate<-ifelse(Nbt==0,NA,-log(Nbt/(Nbt+NFt))/intlen)
	qRate<-ifelse(Nbt==0,NA,-log(Nbt/(Nbt+NbL))/intlen)
	if(plot){
		int.start<-intMat[,1];int.end<-intMat[,2]
		times1<-c(int.start,(int.end+((int.start-int.end)/100)))
		p1<-c(pRate,pRate)[order(times1)]
		q1<-c(qRate,qRate)[order(times1)]
		times1<-sort(times1)
		if(plotLogRates){
			p1[!(p1>0)]<-NA
			q1[!(q1>0)]<-NA
			ylims<-c(min(c(p1,q1),na.rm=TRUE),(max(c(p1,q1),na.rm=TRUE))*1.5)
			plot(times1,p1,type="l",log="y",col=4,lwd=2,lty=5,
				xlim=c(max(times1),max(0,min(times1))),ylim=ylims,
				xlab="Time (Before Present)",ylab="Instantaneous Per-Capita Rate (per Ltu)")
			lines(times1,q1,col=2,lwd=2,lty=2)
			if(showLegend){legend(x="topleft",legend=c("Origination","Extinction"),lty=c(5,2),lwd=2,col=c(4,2))}
		}else{
			ylims<-c(min(c(p1,q1),na.rm=TRUE),(max(c(p1,q1),na.rm=TRUE))*1.2)
			plot(times1,p1,type="l",col=4,lwd=2,lty=5,
				xlim=c(max(times1),max(0,min(times1))),ylim=ylims,
				xlab="Time (Before Present)",ylab="Instantaneous Per-Capita Rate (per Ltu)")
			lines(times1,q1,col=2,lwd=2,lty=2)		
			if(showLegend){legend(x="topleft",legend=c("Origination","Extinction"),lty=c(5,2),lwd=2,col=c(4,2))}
			}
		}
	res<-cbind(intMat,intlen,Nbt,NbL,NFt,NFL,N,pRate,qRate)
	return(invisible(res))
	}

