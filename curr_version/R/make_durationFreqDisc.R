make_durationFreqDisc<-function(timeData,groups=NULL,dropModern=TRUE){
	#this is the multi-parameter maximum likelihood analysis of binned timeData
		#uses a set of binned-interval timeData (just the by-species first and last intervals matrix) 
			#to fit models of different samp probs and ext rates
			#output from binTimeData() can be input directly (only looks at second matrix)
	#NOTE#####
		#UNLIKE getSampProbDisc and getSampRateCont, there are no moving time-windows
		#in fact, this was probably a bad idea to begin with
	#dropModern drops ALL taxa that survive to the modern (i.e. truncated ranges)
	if(length(timeData)==2){	#if a timeList matrix...
		timeList<-timeData
		modernTest<-apply(timeList[[1]],1,function(x) all(x==0))
		if(any(modernTest)){	#if modern present
			if(sum(modernTest)>1){stop("More than one modern interval in timeData??!")}
			#modify the taxon occurrence matrix
			modInt<-which(modernTest)
			newInt<-which(apply(timeList[[1]],1,function(x) x[1]!=0 & x[2]==0))
			if(length(newInt)>1){stop("More than one interval stretching to the modern in timeData??!")}
			if(dropModern){
				modDroppers<-apply(timeList[[2]],1,function(x) x[1]==modInt)
				timeList[[2]]<-timeList[[2]][-modDroppers,]
				if(!is.null(groups)){
					if(dropModern & any(modernTest)){groups<-groups[-modDroppers,]}
					if(nrow(timeList[[2]])!=nrow(groups)){
						stop(paste("number of rows in groups isn't equal to number of taxa in timeData",
							if(dropModern){"after modern taxa are dropped"}))}
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