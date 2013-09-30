multiDiv<-function(data,int.length=1,plot=TRUE,split.int=TRUE,drop.ZLB=TRUE,drop.cryptic=FALSE,
	extant.adjust=0.01,plotLogRich=FALSE,timelims=NULL,int.times=NULL,
	plotMultCurves=FALSE,multRainbow=TRUE,divPalette=NULL){
	#lines up a bunch of taxic or phylo objects and calculates diversity curves simulataneously
		#across all their objects; intuits the type of object without being told
		#it also calculates a "average" median curve and 95% quantile intervals
	#input is a list of dicrete interval or continuous time taxic data or a timetree
		#as in the respective functions
	#output is a list with third objects
		#the first object is a 2-column matrix with interval starts and ends
		#the second object is a matrix 
			#with the measured diversity for all the objects as columns, intervals as rows
	#3rd object consists of a 3col matrix of information related to median curve
		#first column is a per-interval median of the combined diversity curves
		#second and third columns are 95% quantile intervals on that median
	#int.length=1;plot=TRUE;split.int=TRUE;drop.ZLB=TRUE;drop.cryptic=FALSE;plotLogRich=FALSE;timeLims=NULL
	#plotMultCurves=FALSE;multRainbow=TRUE;divPalette=NULL
	#require(ape)
	dclass<-sapply(data,class)	#data classes
	dclass1<-numeric(length(dclass));dclass1[dclass=="matrix"]<-1;
		dclass1[dclass=="list"]<-2;dclass1[dclass=="phylo"]<-3
	if(any(dclass1==0)){stop("Error: Data of Unknown Type")}
	if(is.null(int.times)){
		tblen<-int.length
		#get max and min times for each type
		if(any(dclass1==1)){
			lims1<-sapply(data[dclass1==1],function(x) 
				if(ncol(x)==6){
					c(min(x[,3:4],na.rm=T),max(x[,3:4],na.rm=T))	
				}else{
					c(min(x,na.rm=TRUE),max(x,na.rm=TRUE))}
				)
		}else{lims1<-NA}
		if(any(dclass1==2)){	
			for(i in which(dclass1==2)){
				data[[i]][[1]][data[[i]][[1]][,1]==0,1]<-extant.adjust
				}
			lims2<-sapply(data[dclass1==2],function(x) 
				c(min(x[[1]][max(x[[2]]),]),max(x[[1]][min(x[[2]]),])))
		}else{lims2<-NA}
		if(any(dclass1==3)){
			lims3<-numeric()
			for(i in which(dclass1==3)){
				ttree<-data[[i]]
				if(is.null(ttree$root.time)){
					ntime<-dist.nodes(ttree)[,Ntip(ttree)+1]
					ntime<-max(ntime)-ntime
				}else{
					ntime<-dist.nodes(ttree)[,Ntip(ttree)+1]
					ntime<-ttree$root.time-ntime
					}
				lims3<-c(lims3,c(min(ntime),max(ntime)))
				}
		}else{lims3<-NA}
		end<-min(c(lims1,lims2,lims3),na.rm=TRUE)
		start<-max(c(lims1,lims2,lims3),na.rm=TRUE)
		midtimes<-seq(start+2*tblen,end-2*tblen,by=-tblen)
		midtimes<-midtimes[midtimes>0]
		int.start<-midtimes+(tblen/2);int.end<-midtimes-(tblen/2)
		int.times<-cbind(int.start,int.end)
		if(split.int & any(dclass1==2)){
			#for every single discrete time dataset, bins must be split at their boundaries
			splinters<-sapply(data[dclasss=2],function(x) x[[1]][,1])
			mustSplit<-apply(int.times,1,function(x) any(sapply(splinters,function(y) x[1]>y & x[2]<y)))
			if(any(mustSplit)){
				for(i in which(mustSplit)){
						splitter<-splinters[sapply(splinters[,1],function(y) int.times[i,1]>y & int.times[i,2]<y),1]
						#if(length(splitter)>1){stop("Error: Splitter returning more than one value?!")}
						splitter<-c(int.times[i,1],splitter,int.times[i,2])
						int.times<-rbind(int.times,cbind(splitter[1:(length(splitter)-1)],splitter[2:length(splitter)]))
					}
				int.times<-int.times[-which(mustSplit),]
				int.times<-int.times[order(-int.times[,1]),]
				}
			midtimes<-(int.start+int.end)/2		
			}
		}			
	div<-matrix(,nrow(int.times),1)
	for(i in 1:length(data)){
		if((dclass1[i]==1)){
			divs1<-taxicDivCont(timeData=data[[i]],int.times=int.times,plot=FALSE,drop.cryptic=drop.cryptic)[,3]
			div<-cbind(div,divs1)
			}
		if((dclass1[i]==2)){
			divs2<-taxicDivDisc(timeList=data[[i]],int.times=int.times,split.int=FALSE,plot=FALSE)[,3]
			div<-cbind(div,divs2)
			}
		if((dclass1[i]==3)){
			divs3<-phyloDiv(tree=data[[i]],int.times=int.times,plot=FALSE,drop.ZLB=drop.ZLB)[,3]
			div<-cbind(div,divs3)
			}
		}
	div<-div[,-1]
	colnames(div)<-paste("dataset",1:ncol(div),sep="") 
	#get median curve
	median<-apply(div,1,median)
	q1<-apply(div,1,quantile,probs=0.025)	#the low quantile
	q2<-apply(div,1,quantile,probs=0.975)	#the high quantile
	median.curve<-cbind(median=median,low.95.quantile=q1,high.95.quantile=q2)
	res<-list(int.times=int.times,div.counts=div,median.curve=median.curve)
	if(plot){plotMultiDiv(res,plotLogRich=plotLogRich,timelims=timelims,
		plotMultCurves=plotMultCurves,multRainbow=multRainbow,divPalette=divPalette)}
	return(invisible(res))
	}