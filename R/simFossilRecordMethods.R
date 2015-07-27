#' Methods for Editting or Converting Output from simFossilRecord
#'
#'


#' @name simFossilRecordMethods

#' @details

#' @inheritParams

#' @param

#' @return

 #taxon.id ancestor.id orig.time ext.time still.alive looks.like

#' @aliases

#' @seealso

#' @author 
#' David W. Bapst, inspired by code written by Peter Smits.

#' @references

#' @examples




#' @export
timeSliceFossilRecord<-function(fossilRecord,sliceTime,modern.samp.prob=1){
	#take a fossilRecord data object and cut it at some specific date
	#
	#drop all taxa that originate after the sliceTime
	droppers<-which(sapply(fossilRecord,function(x) x[[1]][3]<sliceTime))
	fossilRecord<-fossilRecord[-droppers]
	#
	#turn all taxa that went extinct after sliceTime so they are still alive
	stillAlive<-which(sapply(fossilRecord,function(x) x[[1]][4]<sliceTime))
	fossilRecord[stillAlive][[1]][4:5]<-c(NA,1)
	#
	#remove all sampling events after sliceTime
	#adjust all dates so cutdate becomes 0
	for(i in 1:length(fossilRecord)){
		#remove all sampling events after sliceTime
		fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]][fossilRecord[[i]][[2]]>=sliceTime]
		#adjust all dates so cutdate becomes 0
		fossilRecord[[i]][[1]][3:4]<-fossilRecord[[i]][[1]][3:4]-sliceTime
		fossilRecord[[i]][[2]]<-fossilRecord[[i]][[2]]-sliceTime
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
	
fossilRecord2fossilTaxa<-function(fossilRecord){
	# a function that transforms a simfossilrecord to a taxa object
	taxaConvert<-t(sapply(fossilRecord,function(x) x[[1]]))	
	rownames(taxaConvert)<-names(fossilRecord)
	return(taxaConvert)
	}

fossilRecord2fossilRanges<-function(fossilRecord, merge.cryptic=TRUE, ranges.only = TRUE){
	# a function that transforms a simfossilrecord to a set of ranges (like from sampleRanges)
		# merge.cryptic = TRUE or FALSE
		# ranges.only or sampling times?
	sampData<-lapply(fossilRecord,function(x) x[[2]]) 
	#merge cryptic taxa
	if(merge.cryptic=TRUE){
		cryptIDs<-sapply(fossilRecord,function(x) x[[1]][6])
		
		}
	rangeData	

	return(rangeData)
	}

# a function that wraps taxa2phylo for simFossilRecord, providing time-scaled tree of sampled taxa
	# merge.cryptic = TRUE or FALSE
	#ala simPaleoTrees:
		# tree<-taxa2phylo(taxa,obs_time=ranges1[,2],plot=plot)	
