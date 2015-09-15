#' Methods for Editing or Converting Output from simFossilRecord
#'
#' These are a set of functions available for manipulating, translating
#' and editing the output data objects from function \code{simFossilRecord}.

#' @name simFossilRecordMethods

#' @details
#' These functions exist to manipulate output from \code{simFossilRecord},
#' particularly so that they can be interfaced with functions in library
#' \code{paleotree} in the same way that output from the
#' (legacy) function \code{simFossilTaxa} was used.

#' @inheritParams sampleRanges

#' @param fossilRecord A list object output by \code{simFossilRecord}, often composed
#' of multiple elements, each of which is data for 'one taxon', with the first
#' element being a distinctive six-element vector composed of numbers, corresponding
#' to the six numbers in a \code{simFossilTaxa} matrix.

#' @param sliceTime The date to slice the \code{simFossilRecord} output at, given
#' in time-units before the modern, on the same scale as the input \code{fossilRecord}.

#' @param shiftRoot4TimeSlice Should the dating of events be shifted, so that the
#' date given for \code{sliceTime} is now 0, or should the dates not be shifted,
#' so that they remain on the same scale as the input? This argument accepts a
#' logical TRUE or FALSE, but also accepts the string \code{"withExtantOnly"},
#' which will only 'shift' the time-scale if living taxa are present, as
#' determined by having ranges that overlap within \code{tolerance} of \code{sliceTime}.

#' @param tolerance A small number which sets a range around the \code{sliceTime} within
#' which taxa will be considered extant.

#' @param modern.samp.prob The probability that a taxon is sampled at the modern time
#' (or, for \code{timeSliceFossilRecord}, the time at which the simulation data is
#' slice). Must be a number between 0 and 1. If 1, all taxa that survive to the modern
#' day (to the \code{sliceTime}) are sampled, if 0, none are.

#' @return
#' Depends on the function.

#' @aliases timeSliceFossilRecord fossilRecord2fossilTaxa fossilRecord2fossilRanges

#' @seealso
#' \code{\link{simFossilRecord}}

#' @author 
#' David W. Bapst.

#' @examples
#'
#' set.seed(444)
#' record <- simFossilRecord(p=0.1, q=0.1, r=0.1, nruns=1,
#' 	nTotalTaxa=c(20,30) ,nExtant=0, plot=TRUE)
#'




#' @rdname simFossilRecordMethods
#' @export
timeSliceFossilRecord<-function(fossilRecord, sliceTime, shiftRoot4TimeSlice=FALSE,
		modern.samp.prob=1, tolerance=10^-4){
	#take a fossilRecord data object and cut it at some specific date
	#
	# CHECKS
	checkResult<-checkFossilRecord(fossilRecord)
	#check shiftRoot4TimeSlice
	shiftPar<-c(TRUE,FALSE,"withExtantOnly")
	shiftRoot4TimeSlice<-shiftPar[pmatch(shiftRoot4TimeSlice,shiftPar)]
	if(is.na(shiftRoot4TimeSlice)){
		stop("shiftRoot4TimeSlice must be a logical or the string 'withExtantOnly'")}
	#
	#drop all taxa that originate after the sliceTime
	droppers<-sapply(fossilRecord,function(x) x[[1]][3]<sliceTime)
	fossilRecord<-fossilRecord[!droppers]
	#
	#remove all sampling events after sliceTime
		for(i in 1:length(fossilRecord)){
			#remove all sampling events after sliceTime
			fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]][fossilRecord[[i]][[2]]>=sliceTime]
		}
	# adjusting time, making taxa extant
		# need to first test if there are extant taxa or not
	isAlive<-sapply(fossilRecord,function(x){
		if(is.na(x[[1]][4])){
			TRUE
		}else{
			(sliceTime-x[[1]][4])>tolerance
		}})
	#browser()
	if(shiftRoot4TimeSlice=="withExtantOnly"){
		if(any(isAlive)){
			shiftRoot4TimeSlice<-TRUE
		}else{
			shiftRoot4TimeSlice<-FALSE
			}
		}
	#
	# if shiftRoot4TimeSlice, then the whole thing shifts so time=0 is slice time
	if(shiftRoot4TimeSlice){	
		#adjust all dates so cutdate becomes 0
		for(i in 1:length(fossilRecord)){
			#adjust all dates so cutdate becomes 0
				#if stillAlive, replace 4:5 with 0,1
			if(isAlive[i]){
				#turn all taxa that went extinct after sliceTime so they are still alive
				fossilRecord[[i]][[1]][3]<-fossilRecord[[i]][[1]][3]-sliceTime
				fossilRecord[[i]][[1]][4:5]<-c(0,1)
			}else{
				fossilRecord[[i]][[1]][3:4]<-fossilRecord[[i]][[1]][3:4]-sliceTime
				}
			fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]]-sliceTime
			}
	# if shiftRoot4TimeSlice=FALSE, then simply replace all extant taxa with
		# LADs at sliceTime and score as extant
	}else{
		for(i in 1:length(fossilRecord)){
			if(isAlive[i]){
				fossilRecord[[i]][[1]][4:5]<-c(sliceTime,1)
				}
			}
		}
	#
	# sample at modern based on modern.samp.prob
	whichExtant<-which(sapply(fossilRecord,function(x) x[[1]][5]==1))
	nLive<-length(whichExtant)
	liveSampled<-as.logical(rbinom(n=nLive, size=1, prob=modern.samp.prob))
	whichSampled<-whichExtant[liveSampled]
	#add sampling event at modern
	for(i in whichSampled){
		fossilRecord[[i]][[2]]<-c(fossilRecord[[i]][[2]],0)
		}
	#
	return(fossilRecord)
	}

#' @rdname simFossilRecordMethods
#' @export
fossilRecord2fossilTaxa<-function(fossilRecord){
	# CHECKS
	checkResult<-checkFossilRecord(fossilRecord)
	# a function that transforms a simfossilrecord to a taxa object
	taxaConvert<-t(sapply(fossilRecord,function(x) x[[1]]))	
	rownames(taxaConvert)<-names(fossilRecord)
	return(taxaConvert)
	}

#' @rdname simFossilRecordMethods
#' @export	
fossilRecord2fossilRanges<-function(fossilRecord, merge.cryptic=TRUE, ranges.only = TRUE){
	# a function that transforms a simfossilrecord to a set of ranges (like from sampleRanges)
		# merge.cryptic = TRUE or FALSE
		# ranges.only or sampling times?
	# CHECKS
	checkResult<-checkFossilRecord(fossilRecord)
	#
	sampData<-lapply(fossilRecord,function(x) x[[2]]) 
	#get sampOcc : separate out the sampling events
	sampOcc<-sapply(fossilRecord,function(x) x[[2]])
	names(sampOcc)<-names(fossilRecord)
	#merge cryptic taxa
	if(merge.cryptic){
		taxonIDs<-sapply(fossilRecord,function(x) x[[1]][1])
		cryptIDs<-sapply(fossilRecord,function(x) x[[1]][6])
		for(i in 1:length(fossilRecord)){
			if(taxonIDs[i]==cryptIDs[i]){
				#if its the original taxon, collect all sampling events
					# for this cryptic complex into one pool
				sampOcc[[i]]<-unlist(c(sampOcc[taxonIDs[i]==cryptIDs]))
				#check that its a vector
				if(is.list(sampOcc[[i]])){
					stop("sampling data for taxa is not coercing correctly to a vector")}
			}else{
				#if its a cryptic taxon that didn't found the complex, erase its data
				sampOcc<-NA
				}
			}
		}
	sampOcc[sapply(sampOcc,length)==0]<-NA
	#convert sampling events to FADs and LADs
	if(ranges.only){
		ranges<-cbind(sapply(sampOcc,max),sapply(sampOcc,min))
		rownames(ranges)<-names(sampOcc)
		colnames(ranges)<-c("FAD","LAD")
		result<-ranges
	}else{
		result<-sampOcc
		}
	return(result)
	}

# don't export
checkFossilRecord<-function(fossilRecord){
	if(any(sapply(fossilRecord,length)!=2)){
		stop("fossilRecord object has taxon entries with more or less than two elements")}
	if(any(sapply(fossilRecord,function(x) length(x[[1]]))!=6)){
		stop("fossilRecord object has taxon entries with more or less than six elements in first element")}		
	return(TRUE)
	}

