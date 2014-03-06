

seqTimeList<-function(timeList,nruns=100,weightSampling=FALSE){
	timeList[[1]]<-as.matrix(timeList[[1]])
	timeList[[2]]<-as.matrix(timeList[[2]])
	#let's try to get sampling rates out of a timeList with overlapping intervals
	#test for overlap
	lap<-any(!apply(timeList[[1]],1,function(x)
		2>sum(x[1]>(timeList[[1]][,2]) & x[2]<(timeList[[1]][,1]))))
	if(!lap){stop("Your timeList doesn't have any overlap??")}
	#let's make a list of valid intervals
		#first, all non-overlapping intervals
	noLap<-which(apply(timeList[[1]],1,function(x) 
		2>sum(x[1]>(timeList[[1]][,2]) & x[2]<(timeList[[1]][,1]))))
	#now let's search iteratively for interval solutions
	nInts<-nTaxa<-numeric()
	savedLists<-list(1:nruns)
	for(i in 1:nruns){
		#set it up!
		valid<-noLap
		invalid<-numeric()
		if(length(valid)>0){
			candidate<-(1:nrow(timeList[[1]]))[-valid]
			}else{
			candidate<-(1:nrow(timeList[[1]]))
			}
		#LOOP
		while(length(candidate)>1){
			#next, let's pull out a random candidate interval
			if(weightSampling){
				weights<-order(apply(timeList[[1]][candidate,,drop=FALSE],1,diff))
				weights<-weights/sum(weights)
				draw<-sample(1:length(candidate),1,prob=weights)
			}else{
				draw<-sample(1:length(candidate),1)
				}
			#remove from candidate, add to valid
			valid<-c(valid,candidate[draw])		
			drawInt<-timeList[[1]][candidate[draw],]
			candidate<-candidate[-draw]
			#add all overlapping to invalid and remove from candidate
			candInt<-timeList[[1]][candidate,,drop=FALSE]
			badInt<-apply(candInt,1,function(x) drawInt[1]>x[2] & drawInt[2]<x[1])	
			invalid<-c(invalid,candidate[badInt])
			candidate<-candidate[!badInt]
			#
			#newInts<-timeList[[1]][valid,]
			#if(any(!apply(newInts,1,function(x)
			#	2>sum(x[1]>(newInts[,2]) & x[2]<(newInts[,1]))))){
			#		stop("WHAT!")}
			}
		if(length(candidate)==1){valid<-c(valid,candidate)}
		valid<-sort(valid)
		newInts<-timeList[[1]][valid,]
		# test that there be no overlaps
		if(any(!apply(newInts,1,function(x)	
			2>sum(x[1]>(newInts[,2]) & x[2]<(newInts[,1]))))){stop("WHAT!")}
		#
		validTaxa<-apply(timeList[[2]],1,function(x) any(x[1]==valid) & any(x[2]==valid))
		cropTaxa<-timeList[[2]][validTaxa,]
		newTaxa<-cbind(sapply(cropTaxa[,1],function(x) which(x==valid)),
		sapply(cropTaxa[,2],function(x) which(x==valid)))
		nInts[i]<-nrow(newInts)
		nTaxa[i]<-nrow(newTaxa)
		savedLists[[i]]<-list(intTimes=newInts,taxonTimes=newTaxa)
		}
	results<-list(nIntervals=nInts,nTaxa=nTaxa,timeLists=savedLists)
	return(results)
	}