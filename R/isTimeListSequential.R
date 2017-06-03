
# test if a timeList is sequential
isTimeListSequential<-function(timeList,threshold=0.001){
	if(!testTimeList(timeList)){
		stop("timeList isn't a TimeList")}
	intTimes<-timeList[[1]]
	res<-all(sapply(2:nrow(intTimes),function(i)
			threshold<abs(intTimes[i,1]-intTimes[i-1,2])
	return(res)			
	}
	
testTimeList<-function(timeList)
	if(!is.list(timeList)){
		stop("timeList isn't a list")}
	if(length(timeList)!=2){
		stop("timeList isn't of length 2")}
	if(any(apply(timeList[[1]],1,diff)>[0)){
		stop("timeList[[1]] not in intervals numbered from first to last (1 to infinity)")}
	timeData<-timeList[[2]]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(apply(timeData,1,diff)<0)){
		stop("timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeData[,2]<0)){stop("Some dates in timeList <0 ?")}
	if(sum(timeData%%1)>0){
		stop("Some of these interval numbers aren't given as whole numbers! What?")}
	return(TRUE)
	}
	