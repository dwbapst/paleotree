#' Models of Sampling and Extinction for Taxonomic Duration Datasets
#'
#' These functions construct likelihood models of the observed frequency
#' of taxon durations, given either in discrete (make_durationFreqDisc) 
#' or continuous time (make_durationFreqCont). These models can then be
#' constrained using functions available in this package and/ord analyzed 
#' with commonly used optimizing functions.
#'
#' @details
#' These functions effectively replace two older functions in paleotree,
#' \code{\link{getSampRateCont}} and \code{\link{getSampRateDisc}}, which
#' are otherwise retained in paleotree for historical purposes. The
#' functions here do not offer the floating time interval options of
#' their older siblings, but do allow for greater flexibility in defining
#' constrains on parameter values. Differences in time intervals, or any
#' other conceivable discrete differences in parameters, can be modeled
#' using the generic \code{groups} argument in these functions.
#'
#' These functions use likelihood functions presented by Foote (1997).
#' These analyses are ideally applied to data from single stratigraphic section
#' but potentially are applicable to regional or global datasets, although the
#' behavior of those datasets is less well understood.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero and older dates are 'larger'. On
#' the contrary, relative time is in intervals with non-zero integers that
#' increase sequentially beginning with 1, from earliest to oldest.
#'
#' For make_durationFreqDisc, the intervals in timeList should be
#' non-overlapping sequential intervals of roughly equal length. These
#' should be in relative time as described above, so the earliest interval
#" should be 1 and the numbering should increase as the intervals go up with
#' age. If both previous statements are true, then differences in interval
#' numbers will represent the same rough difference in the absolute timing
#" of those intervals. For example, a dataset where all taxa are listed from
#' a set of sequential intervals of similar length, such as North American
#' Mammal assemblage zones, microfossil faunal zones or graptolite biozones
#" can be given as long as they are correctly numbered in sequential order
#" in the input. As a counter example, a dataset which includes taxa resolved
#" only to intervals as wide as the whole Jurassic and taxa resolved to
#" biozones within the Jurassic should not be included in the same input.
#' Drop taxa from less poorly resolved intervals from such datasets if you
#" want to apply this function, as long as this retains a large enough sample
#" of taxa listed from the sequential set of intervals. 
#' 
#' Please check that the optimizer function you select actually converges. The
#' likelihood surface can be very flat in some cases, particularly for small
#' datasets (<100 taxa). If the optimizer does not converge, consider
#' increasing iterations or changing the starting values.
#'

#' @param timeData Two-column matrix of per-taxon first and last occurrence
#' given in continuous time, relative to the modern (i.e. older dates are also
#' the 'larger' dates).

#' @param timeList A 2 column matrix with the first and last occurrences of taxa
#' given in relative time intervals (i.e. ordered from first to last). If a list
#' of length two is given for timeData, such as would be expected if the output 
#' of binTimeData was directly input, the second element is used. See details.

#' @param groups Either NULL (the default) or matrix with the number of rows equal
#' to the number of taxa and the number of columns equal to the number of 'systems'
#' of categories for taxa. Taxonomic membership in different groups is indicated
#' by numeric values. For example, a dataset could have a 'groups' matrix composed
#' of a column representing thin and thick shelled taxa, coded 1 and 2 respectively,
#' while the second column indicates whether taxa live in coastal, shelfal or deep
#' marine settings, coded 1-3 respectively. Different combinations of groups will
#' be treated as having independent sampling and extinction parameters in the
#' default analysis, for example, thinly-shelled deep marine species will have
#' separate parameters from thinly-shelled coastal species. Grouping systems could
#' also represent temporal heterogeneity, for example, categorizing Paleozoic versus
#' Mesozoic taxa. If groups are NULL (the default), all taxa are assumed to be of
#' the same group with the same parameters.

#' @param drop.extant Drops all extant taxa from a dataset.

#' @param threshold The smallest allowable duration (i.e. the measured difference in
#' the first and last occurrence dates for a given taxon). Durations below this size 
#' will be treated as "one-hit" sampling events.

#' @return 
#' A function of class "paleotreeFunc", which takes vector equal to the number
#' of parameters and returns the *negative* log likelihood (for use with optim and
#' similar optimizing functions, which attempt to minimize support values). See the
#' functions listed at 

#' @aliases make_durationFreqCont make_durationFreqDisc

#' @author David W. Bapst

#' @seealso
#' See the original implementation of these methods at 
#' \code{\link{getSampRateCont}} and \code{\link{getSampProbDisc}}. 
#'
#' Also see \code{\link{freqRat}}, \code{\link{sRate2sProb}},
#' \code{\link{qsRate2Comp}} \code{\link{sProb2sRate}} and \code{\link{qsProb2Comp}}.

#' @references 
#' Foote, M. 1997 Estimating Taxonomic Durations and Preservation
#' Probability. \emph{Paleobiology} \bold{23}(3):278--300.
#' 
#' Foote, M., and D. M. Raup. 1996 Fossil preservation and the stratigraphic
#' ranges of taxa. \emph{Paleobiology} \bold{22}(2):121--140.

#' @examples
#' 

#' @name durationFreq
#' @rdname durationFreq
#' @export
make_durationFreqCont<-function(timeData,groups=NULL,drop.extant=TRUE,threshold=0.01){
	#this is the multi-parameter maximum likelihood analysis of binned timeData
		#uses a set of binned-interval timeData (just the by-species first and last intervals matrix) 
			#to fit models of different samp probs and ext rates
			#output from binTimeData() can be input directly (only looks at second matrix)
	#NOTE####
		#UNLIKE getSampProbDisc and getSampRateCont, there are no moving time-windows
		#in fact, this was probably a bad idea to begin with
	#drop.extant drops ALL taxa that survive to the modern (i.e. truncated ranges)
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	timeData<-timeData[!is.na(timeData[,1]) & !is.na(timeData[,2]),]
	modernTest<-apply(timeData,1,function(x) all(x==0))
	if(any(modernTest)){	#if modern present
		#modify the taxon occurrence matrix
		if(drop.extant){
			modDroppers<-timeData[,2]==0
			timeData<-timeData[-modDroppers,]
			if(!is.null(groups)){
				if(drop.extant & any(modernTest)){groups<-groups[-modDroppers,]}
				if(nrow(timeData)!=nrow(groups)){
					stop(paste("number of rows in groups isn't equal to number of taxa in timeData",
						if(drop.extant){"after modern taxa are dropped"}))}
				}
			}
		}
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("Error: timeData is not in time relative to modern (decreasing to present)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData <0 ?")}
	#get the dataset
	dur<-(timeData[,1]-timeData[,2])
	#THRESHOLD DETERMINES RANGES TOO SMALL TO BE CONSIDERED NOT ONE-TIMERS
	dur[dur<threshold]<-0
	#define parnames
	parnames<-c("q","r")
	if(is.null(groups)){groups2<-matrix(1,length(dur),1)}else{groups2<-groups}
	for(i in 1:ncol(groups2)){
		parnames<-as.vector(sapply(parnames,function(x) paste(x,unique(groups2[,i]),sep=".")))
		}
	groupings<-unique(groups2)
	ngroup<-nrow(groupings)
	#break parnames into a character matrix
	breakNames<-t(sapply(parnames,function(x) unlist(strsplit(x,split=".",fixed=TRUE))))
	#NEED TO FIGURE OUT parbounds
	lowerBound<-rep(0.001,length(parnames))
	upperBound<-rep(5,length(parnames))
	parbounds<-list(lowerBound,upperBound)
	#
	logL<-function(par){
		if(length(par)!=length(parnames)){stop("Number of input parameters is not equal to number of parnames")}
		logLsum<-numeric()
		for(i in 1:ngroup){
			#if(is.null(groups)){selector<-rep(TRUE,length(par))
			#	}else{selector<-apply(breakNames,1,function(x) x[-(1:2)]==groupings[i,])}
			selector<-apply(breakNames,1,function(x) x[-1]==groupings[i,])
			parnew<-par[selector]
			breaknew<-breakNames[selector,]
			q<-parnew[breaknew[,1]=="q"]
			r<-parnew[breaknew[,1]=="r"]
			dur1<-dur[selector]
			ft<-ifelse(dur1==0,log(q/(r+q)),log(q*r*exp(-q*dur1)/(r+q)))
			logLsum[i]<-(sum(ft))
			}		
		res<-(-sum(logLsum))
		return(unname(res))
		}
	#make into a paletree likelihood function
	logL<-make_paleotreeFunc(logL,parnames,parbounds)
	return(logL)
	}

#' @rdname durationFreq
#' @export	
make_durationFreqDisc<-function(timeList,groups=NULL,drop.extant=TRUE){
	#this is the multi-parameter maximum likelihood analysis of binned timeData
		#uses a set of binned-interval timeData (just the by-species first and last intervals matrix) 
			#to fit models of different samp probs and ext rates
			#output from binTimeData() can be input directly (only looks at second matrix)
	#NOTE#####
		#UNLIKE getSampProbDisc and getSampRateCont, there are no moving time-windows
		#in fact, this was probably a bad idea to begin with
	#drop.extant drops ALL taxa that survive to the modern (i.e. truncated ranges)
	timeData<-timeList	#because i'm an idiot and i dislike having to rename arguments?
	if(length(timeData)==2){	#if a timeList matrix...
		timeList<-timeData
		modernTest<-apply(timeList[[1]],1,function(x) all(x==0))
		if(any(modernTest)){	#if modern present
			if(sum(modernTest)>1){stop("More than one modern interval in timeData??!")}
			#modify the taxon occurrence matrix
			modInt<-which(modernTest)
			newInt<-which(apply(timeList[[1]],1,function(x) x[1]!=0 & x[2]==0))
			if(length(newInt)>1){stop("More than one interval stretching to the modern in timeData??!")}
			if(drop.extant){
				modDroppers<-apply(timeList[[2]],1,function(x) x[1]==modInt)
				timeList[[2]]<-timeList[[2]][-modDroppers,]
				if(!is.null(groups)){
					if(drop.extant & any(modernTest)){groups<-groups[-modDroppers,]}
					if(nrow(timeList[[2]])!=nrow(groups)){
						stop(paste("number of rows in groups isn't equal to number of taxa in timeData",
							if(drop.extant){"after modern taxa are dropped"}))}
					}
				}
			#change all modInt references to the prior int in the taxon appearance matrix
			timeList[[2]]<-apply(timeList[[2]],2,sapply,function(x) if(x==modInt){newInt}else{x})
			}
		timeData<-timeList[[2]]
		}
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(timeData,1,diff)<0)){
		stop("Error: timeData / timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Error: Some dates in timeData / timeList[[2]] <0 ?")}
	if(sum(timeData%%1)>0){stop("Error: Some of these interval numbers aren't given as whole numbers! What?")}
	#get the dataset
	dur<-apply(timeData,1,diff)+1
	#timeData1<-max(timeData)-timeData+1
	#FO<-timeData1[,1];LO<-timeData1[,2]	#not needed
	#define parnames
	parnames<-c("q","R")
	#BAD
	#if(!is.null(groups)){
	#	for(i in 1:ncol(groups)){
	#		parnames<-as.vector(sapply(parnames,function(x) paste(x,unique(groups[,i]),sep=".")))
	#		}
	#	groupings<-unique(groups)
	#	}
	#ngroup<-ifelse(is.null(groups),1,nrow(groupings))
	#NEW
	if(is.null(groups)){groups2<-matrix(1,length(dur),1)}else{groups2<-groups}
	for(i in 1:ncol(groups2)){
		parnames<-as.vector(sapply(parnames,function(x) paste(x,unique(groups2[,i]),sep=".")))
		}
	groupings<-unique(groups2)
	ngroup<-nrow(groupings)
	#
	#break parnames into a character matrix
	breakNames<-t(sapply(parnames,function(x) unlist(strsplit(x,split=".",fixed=TRUE))))
	#NEED TO FIGURE OUT parbounds
	lowerBound<-rep(0.001,length(parnames))
	upperBound<-rep(100,length(parnames))
	upperBound[breakNames[,1]=="R"]<-1
	parbounds<-list(lowerBound,upperBound)
	logL<-function(par){
		if(length(par)!=length(parnames)){stop("Number of input parameters is not equal to number of parnames")}
		logLsum<-numeric()
		for(i in 1:ngroup){
			#if(is.null(groups)){selector<-rep(TRUE,length(par))
			#	}else{selector<-apply(breakNames,1,function(x) x[-(1:2)]==groupings[i,])}
			selector<-apply(breakNames,1,function(x) x[-1]==groupings[i,])
			parnew<-par[selector]
			breaknew<-breakNames[selector,]
			q<-parnew[breaknew[,1]=="q"]
			R<-parnew[breaknew[,1]=="R"]
			dur1<-dur[selector]
			TMax<-max(dur1)
			Ti<-1:TMax	#unique durations
			N<-sapply(Ti,function(t) sum(t==dur1))		#number with those durations
			PDT<-exp(-q*(Ti-1))-exp(-q*Ti)
			Rex<-c(1,rep(2,TMax-1))
			ft<- sapply(Ti,function(t) sum(PDT[t:TMax]*((R^Rex[t])*(Ti[t:TMax]-t+1)*((1-R)^(Ti[t:TMax]-t)))))
			ft<-log(ft/sum(ft))	#divide by sum; same as normalizing by Pp acc. to Foote (really? COOL!)
			logLsum[i]<-sum((N*ft)[!is.infinite(ft)])
			}		
		res<-(-sum(logLsum))
		return(unname(res))
		}
	#make into a paletree likelihood function
	logL<-make_paleotreeFunc(logL,parnames,parbounds)
	return(logL)
	}