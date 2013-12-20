#' @details

#' @inheritParam

#' @param

#' @return

#' @aliases

#' @seealso

#' @references

#' @examples

make_durationFreqCont<-function(timeData,groups=NULL,dropModern=TRUE,threshold=0.01){
	#this is the multi-parameter maximum likelihood analysis of binned timeData
		#uses a set of binned-interval timeData (just the by-species first and last intervals matrix) 
			#to fit models of different samp probs and ext rates
			#output from binTimeData() can be input directly (only looks at second matrix)
	#NOTE####
		#UNLIKE getSampProbDisc and getSampRateCont, there are no moving time-windows
		#in fact, this was probably a bad idea to begin with
	#dropModern drops ALL taxa that survive to the modern (i.e. truncated ranges)
	if(class(timeData)!="matrix"){if(class(timeData)=="data.frame"){timeData<-as.matrix(timeData)
		}else{stop("Error: timeData not of matrix or data.frame format")}}
	timeData<-timeData[!is.na(timeData[,1]) & !is.na(timeData[,2]),]
	modernTest<-apply(timeData,1,function(x) all(x==0))
	if(any(modernTest)){	#if modern present
		#modify the taxon occurrence matrix
		if(dropModern){
			modDroppers<-timeData[,2]==0
			timeData<-timeData[-modDroppers,]
			if(!is.null(groups)){
				if(dropModern & any(modernTest)){groups<-groups[-modDroppers,]}
				if(nrow(timeData)!=nrow(groups)){
					stop(paste("number of rows in groups isn't equal to number of taxa in timeData",
						if(dropModern){"after modern taxa are dropped"}))}
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