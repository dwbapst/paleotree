horizonSampRate<-function(sampOcc=NULL,d=NULL,n=NULL){
	#for estimating sampling rate from continuous time data
		#taken from Solow & Smith, 1997, Paleobiology
	#input is a list, with each element a taxon, consisting of n sampling occurrences
		#just like sampleRanges(data,ranges.only=FALSE)
	if(is.null(d) & is.null(n)){
		if(!is.list(filtered)){
			stop("Error: sampOcc isn't a list of species occurrences")}
		filtered<-sampOcc[!is.na(sampOcc)]
		n<-sapply(filtered,length)
		d<-sapply(filtered,max)-sapply(filtered,min)
		}
	if(is.null(d) | is.null(n)){
		stop("Error: d (duration) and n (# sampling events) have to be both supplied, if one is given")}
	sampRate<-(sum(n-1)^2)/(sum(d)*sum(n))
	names(sampRate)<-"sampRate"
	return(sampRate)
	}