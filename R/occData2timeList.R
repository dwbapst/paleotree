#' Converting Occurrences Data to a timeList Data Object
#'
#' This function converts occurrence data, given as a list where each element
#' is a different taxon's occurrence table (containing minimum and maximum ages
#' for each occurrence), to the 'timeList' format, consisting of a list consisting
#' of a matrix of lower and upper age bounds for intervals and a second matrix
#' recording the interval in which taxa first and last occur in the given dataset.

#the following nonsense is stuff from a blog post that I'm going to turn into documentation

# turn these PBDB occurrences into a `paleotree` `timelist` data format. As you might recall from [some](http://nemagraptus.blogspot.com/2013/06/a-tutorial-to-cal3-time-scaling-using.html) [previous posts](http://nemagraptus.blogspot.com/2015/02/how-do-we-treat-fossil-age-data-dates.html), `timeList` data objects are lists composed of two matrices, the first matrix giving the start and end times of intervals, and the second matrix giving the intervals in which taxa are first and last observed.

#Now, I imagine there might be more reasons to need to convert occurrence data to a `timeList` object than just those involving the PBDB (such as output from simulation function `sampleRanges`), so I'll allow this function to look either for (a) variables named `early_age` and `late_age` as under 'pbdb' vocab, or, (b) if each element of the list is a two-column matrix, use each pair of values as the earliest and latest time bounds for the listed occurrences.

#' @details
#' This function should translate taxon-sorted occurrence data (i.e. datasets output by \code{\link{taxonSortPBDBocc}} 

#' @param occList A list where every element is a table of occurrence data for a different taxon,
#' such as that returned by \code{\link{taxonSortPBDBocc}}. The occurrence data can be either a 
#' two-column matrix composed of the lower and upper age bounds on each taxon occurrence, or has
#' two named variables which match either the 

#' @param intervalType Must be either "dateRange" (the default) or "occRange". If "dateRange",
#' the bounds on the first appearances are the span between the oldest upper and lower bounds
#' of the occurrences, and the bounds on the last appearances are the span between the youngest
#' upper and lower bounds across all occurrences. This is guaranteed to provide the smallest
#' bounds on the first and last appearances, and was suggested by Jon Marcot. The "occRange" option
#' instead returns the smallest bounds among (a) the oldest occurrences for the first appearance
#' (i.e. all occurrences with their lowest bound at the oldest lower age bound), and (b) the
#' youngest occurrences for the last appearance (i.e. all occurrences with their uppermost bound
#' at the youngest upper age bound).

#' @return
#' Returns a standard timeList data object, as used by many other paleotree functions, like \code{\link{bin_timePaleoPhy}}, \code{\link{bin_cal3TimePaleoPhy}} and \code{\link{taxicDivDisc}}

#' @seealso
#' \code{\link{taxonSortPBDBocc}} and the example graptolite dataset at \code{\link{graptPBDB}}

#' @author 
#' David W. Bapst, with the 'dateRange' algorithm suggested by Jon Marcot.

#' @examples
#' data(graptPBDB)
#' 
#' graptOccSpecies<-taxonSortPBDBocc(graptOccPBDB,rank="species",onlyFormal=FALSE)
#' graptTimeSpecies<-occData2timeList(occList=graptOccSpecies)
#' 
#' head(graptTimeSpecies[[1]])
#' head(graptTimeSpecies[[2]])
#' 
#' graptOccGenus<-taxonSortPBDBocc(graptOccPBDB,rank="genus",onlyFormal=FALSE)
#' graptTimeGenus<-occData2timeList(occList=graptOccGenus)
#' 
#' layout(1:2)
#' taxicDivDisc(graptTimeSpecies)
#' taxicDivDisc(graptTimeGenus)
#' 
#' # the default interval calculation is "dateRange"
#' # let's compare to the other option, "occRange"
#' 	#for species
#' 
#' graptOccRange<-occData2timeList(occList=graptOccSpecies, intervalType="occRange")
#' 
#' #we would expect no change in the diversity curve
#' 	#because there are only changes in th
#' 		#earliest bound for the FAD
#' 		#latest bound for the LAD
#' #so if we are depicting ranges within maximal bounds
#' 	#dateRanges has no effect
#' layout(1:2)
#' taxicDivDisc(graptTimeSpecies)
#' taxicDivDisc(graptOccRange)
#' #yep, identical
#' 
#' #so how much uncertainty was gained by using dateRange?
#'
#' # write a simple function for getting uncertainty in first and last
#' 		# appearance dates from a timeList object
#' sumAgeUncert<-function(timeList){
#' 	fourDate<-timeList2fourDate(timeList)
#' 	perOcc<-(fourDate[,1]-fourDate[,2])+(fourDate[,3]-fourDate[,4])
#' 	sum(perOcc)
#' 	}
#'
#' #total amount of uncertainty in occRange dataset
#' sumAgeUncert(graptOccRange)
#' #total amount of uncertainty in dateRange dataset
#' sumAgeUncert(graptTimeSpecies)
#' #the difference
#' sumAgeUncert(graptOccRange)-sumAgeUncert(graptTimeSpecies)
#' #as a proportion
#' 1-(sumAgeUncert(graptTimeSpecies)/sumAgeUncert(graptOccRange))
#' 
#' #a different way of doing it
#' dateChange<-timeList2fourDate(graptTimeSpecies)-timeList2fourDate(graptOccRange)
#' apply(dateChange,2,sum)
#' #total amount of uncertainty removed by dateRange algorithm
#' sum(abs(dateChange))
#' 
#' layout(1)


#' @name occData2timeList
#' @rdname occData2timeList
#' @export
occData2timeList<-function(occList,intervalType="dateRange"){
		#intervalType="dateRange"
	#the following is all original, though inspired by paleobioDB code
	#need checks
	#is occList a list
	if(!is.list(occList)){stop("occList is not a list?")}
	#are all elements of occList matrices
	if(!(all(sapply(occList,is,class2="data.frame")) | all(sapply(occList,is,class2="matrix")))){
		stop("All elements of occList must be all of type data.frame or type matrix")}
	#pull first list entry as an example to check with
	exOcc<-occList[[1]]
	#test if all occList entries have same number of columns
	if(!all(sapply(occList,function(x) ncol(x)==ncol(exOcc)))){
		stop("Not all occList entries have same number of columns?")}
	#will assume all data given in this manner either has columns named "" and ""
		#or has only two columns
	if(any(colnames(exOcc)=="early_age") & any(colnames(exOcc)=="late_age")){
		ageSelector<-c("early_age","late_age")
 	}else{
		if(any(colnames(exOcc)=="eag") & any(colnames(exOcc)=="lag")){
			if(ncol(exOcc)!=2){
				ageSelector<-1:2
			}else{
				stop("Data is not a list of two-column matrics *and* lacks named age columns (from the PBDB)")
				}
		}else{
			ageSelector<-c("early_age","late_age")
			}
		}
	#
	if(!any(intervalType==c("occRange","dateRange"))){stop("intervalType must be one of 'dateRange' or 'occRange'")}
	#get intervals in which taxa appear
	taxaInt<-lapply(occList,function(x) x[,ageSelector])
	if(intervalType=="occRange"){
		#get earliest and latest intervals the taxa appear in
		taxaMin<-lapply(taxaInt,function(x) x[x[,1]==max(x[,1]),,drop=FALSE])
		taxaMax<-lapply(taxaInt,function(x) x[x[,2]==min(x[,2]),,drop=FALSE])
		#get smallest intervals
		taxaMin<-t(sapply(taxaMin,function(x) if(nrow(x)>1){
			x[which((-apply(x,1,diff))==min(-apply(x,1,diff)))[1],]
			}else{x}))
		taxaMax<-t(sapply(taxaMax,function(x) if(nrow(x)>1){
			x[which((-apply(x,1,diff))==min(-apply(x,1,diff)))[1],]
			}else{x}))
		#transpose and remove weird list attribute
		#taxaMin<-t(taxaMin)
		#taxaMax<-t(taxaMax)
		taxaMin<-cbind(unlist(taxaMin[,1]),unlist(taxaMin[,2]))
		taxaMax<-cbind(unlist(taxaMax[,1]),unlist(taxaMax[,2]))
		}
	if(intervalType=="dateRange"){
		#this algorithm suggested by Jon Marcot
			#should be smallest possible interval for FAD and LAD
		#get the earliest early age and earliest latest age as taxaMin (FAD)
		taxaMin<-t(sapply(taxaInt,function(x) apply(x,2,max)))
		#get the earliest early age and earliest latest age as taxaMax (LAD)
		taxaMax<-t(sapply(taxaInt,function(x) apply(x,2,min)))
		}
	fourDate<-cbind(taxaMin,taxaMax)
	timeList<-fourDate2timeList(fourDate)
	rownames(timeList$taxonTimes)<-names(occList)
	return(timeList)
	}
	