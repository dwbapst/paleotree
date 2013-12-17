#' perCapitaRates
#' 
#' Calculates and plots per-capita origination and extinction rates from
#' sequential discrete-time taxon ranges, following Foote (2000).
#' 
#' Background: Diversity curves are plots of species/taxon/lineage richness
#' over time for a particular group of organisms. For paleontological studies,
#' these are generally based on per-taxon range data while more recently in
#' evolutionary biology, molecular phylogenies have been used to calculate
#' lineage-through-time plots (LTTs). Neither of these approaches are without
#' their particular weaknesses; reconstructing the true history of biodiversity
#' is a difficult task no matter what data is available.
#' 
#' The diversity curves produced by these functions will always measure
#' diversity within binned time intervals (and plot them as rectangular bins).
#' For continuous-time data or phylogenies, one could decrease the int.length
#' used to get what is essentially an 'instantaneous' estimate of diversity.
#' This is warned against, however, as most historical diversity data will have
#' some time-averaging or uncertain temporal resolution and thus is probably
#' not finely-resolved enough to calculate instantaneous estimates of
#' diversity.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' As diversity is counted within binned intervals, the true standing diversity
#' may be somewhat lower than the measured (observed) quantity, particularly if
#' intervals are longer than the mean duration of taxa is used. This will be an
#' issue with all diversity curve functions, but particularly the discrete-time
#' variant. For diversity data in particularly large discrete time intervals,
#' plotting this data in smaller bins which do not line up completely with the
#' original intervals will create a 'spiky' diversity curve, as these smaller
#' intersecting bins will have a large number of taxa which may have been
#' present in either of the neighboring intervals. This will give these small
#' bins an apparently high estimated standing diversity. This artifact is
#' avoided with the default setting split.int=TRUE, which will split any input
#' or calculated intervals so that they start and end at the boundaries of the
#' discrete-time range bins.
#' 
#' The timeList object should be a list composed of two matrices, the first
#' matrix giving by-interval start and end times (in absolute time), the secont
#' matrix giving the by-taxon first and last appearances in the intervals
#' defined in the first matrix, numbered as the rows. Absolute time should be
#' decreasing, while the intervals should be numbered so that the number
#' increases with time. Unlike getSampProbDisc, the intervals can be
#' overlapping. See the documentation for binTimeData for more information.
#' 
#' phyloDiv will resolve polytomies to be dichotomous nodes seperated by
#' zero-length branches prior to calculating the diversity curve. There is no
#' option to alter this behavior, but it should not affect the use of the
#' function because the addition of the zero-length branches should produce an
#' identical diveristy history as a polytomy. phyloDiv will also drop
#' zero-length terminal branches, as with the function dropZLB. This the
#' default behavior for the function but can be turned off by setting the
#' argument drop.zlb to FALSE.
#' 
#' @param timeList A list giving interval data. See details below.
#' @param plot If true, diversity curve is plotted
#' @param plotLogRich If true, taxic diversity plotted on log scale
#' @return These functions will invisibly return a three-column matrix, where
#' the first two columns are interval start and end times and the third column
#' is the number of taxa/lineages counted in that interval.
#' @author David W. Bapst
#' @seealso \code{\link{DiversityCurves}}, \code{\link{binTimeData}}
#' 
#' There are several different functions for traditional LTT plots
#' (phylogenetic diversity curves), such as the function
#' ,\code{\link{ltt.plot}} in the package ape, the function \code{ltt} in the
#' package phytools, the function \code{plotLtt} in the package laser and the
#' function \code{LTT.average.root} in the package TreeSim.
#' @examples
#' 
#' #Simulate some fossil ranges with simFossilTaxa
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=80,maxtaxa=100,maxtime=1000,maxExtant=0)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #Now let's use binTimeData() to bin in intervals of 5 time units
#' rangesDisc <- binTimeData(rangesCont,int.length=5)
#' #and get the per-capita rates
#' perCapitaRates(rangesDisc)
#' #on a log scale
#' perCapitaRates(rangesDisc,plotLogRates=TRUE)
#' 
#' #get mean and median per-capita rates
#' res<-perCapitaRates(rangesDisc,plot=FALSE)
#' apply(res[,c("pRate","qRate")],2,mean,na.rm=TRUE)
#' apply(res[,c("pRate","qRate")],2,median,na.rm=TRUE)
#' 
#' #with modern taxa
#' set.seed(444)
#' taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=80,maxtaxa=100,maxtime=1000)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5,,modern.samp.prob=1)
#' #Now let's use binTimeData() to bin in intervals of 5 time units
#' rangesDisc <- binTimeData(rangesCont,int.length=5)
#' #and now get per-capita rates
#' perCapitaRates(rangesDisc)
#' 
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

