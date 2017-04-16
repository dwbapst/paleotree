#' Construct A Block of Tip Age Calibrations for Use with Tip-Dating Analyses in MrBayes
#' 
#' Takes a set of tip ages (in several possible forms, see below),
#' and outputs a set of tip age calibrations
#' for use with tip-dating analyses (sensu Zhang et al., 2016)
#' in the popular phylogenetics program \emph{MrBayes}.
#' These calibrations are printed as a set of character strings, as well as a 
#' line placing an offset exponential prior on the tree age, either
#' printed in the R console or in a named text file, which can be used as
#' commands in the \emph{MrBayes} block of a NEXUS file for use with 
#' (you guessed it!) \emph{MrBayes}.

#' @details
#' Beware: some combinations of arguments might not make sense for your data.
#' 
#' (But that's always true, isn't it?)

#' @param tipTimes This input may be either a timeList object (i.e. a list of length 2, 
#' composed of a table of interval upper and lower time boundaries (i.e. earlier and latter bounds), and 
#' a table of first and last intervals for taxa) or a matrix with rownames
#' for taxa as you want listed in the MrBayes block, with either one, two
#' or four columns containing ages (respectively) for point occurences with
#' precise dates (for a single column), uncertainty bounds on a point occurrence
#' (for two columns), or uncertainty bounds on the first and
#' last occurrence (for four columns). Note that precise first and last occurrence
#' dates should not be entered as a two column matrix, as this will instead be interpreted
#' as uncertainty bounds on a single occurrence. Instead, either select which you want to
#' use for tip-dates and give a 1-column matrix, or repeat (and collate) the columns, so that
#' the first and last appearances has uncertainty bounds of zero.

#' @param whichAppearance Which appearance date of the taxa should be used:
#' their \code{'first'} or their \code{'last'} appearance date? The default
#' option is to use the 'first' appearance date. Note that use of the last
#' appearance date means that tips will be constrained to occur before their
#' last occurrence, and thus could occur long after their first occurrence (!).

#' @param ageCalibrationType This argument decides how age calibrations are defined, 
#' and currently allows for four options: \code{"fixedDateEarlier"} which fixes tip
#' ages at the earlier (lower) bound for the selected age of appearance (see argument
#' \code{whichAppearance} for how that selection is made), \code{"fixedDateLatter"}
#' which fixes the date to the latter (upper) bound of the selected age of appearance,
#' \code{"fixedDateRandom"} which fixes tips to a date that is randomly drawn from a
#' uniform distribution bounded by the upper and lower bounds on the selected age of
#' appearance, or (the recommended option) \code{"uniformRange"} which places a uniform
#' prior on the age of the tip, bounded by the latest and earliest (upper and lower)
#' bounds on the the selected age.

#' @param treeAgeOffset A parameter given by the user controlling the offset 
#' between the minimum and expected tree age prior. mean tree age for the
#' offset exponential prior on tree age will be set to the minimum tree age, 
#' plus this offset value. Thus, an offset of 10 million years would equate to a prior
#' assuming that the expected tree age is around 10 million years before the minimum age.

#' @param minTreeAge if \code{NULL} (the default), then minTreeAge will
#' be set as the oldest date among the tip age used (those used being
#' determine by user choices (or oldest bound on a tip age). Otherwise,
#' the user can supply their own minimum tree, which must be greater than
#' whatever the oldest tip age used is.

#' @param file Filename (possibly with path) as a character string
#' to a file which will be overwritten with the output tip age calibrations.
#' If \code{NULL}, tip calibration commands are output to the console.

#' @return
#' If argument \code{file} is \code{NULL}, then the tip age commands
#' are ouput as a series of character strings.

#' @author
#' David W. Bapst. This code was produced as part of a project 
#' funded by National Science Foundation grant EAR-1147537 to S. J. Carlson.

#' @references
#' Zhang, C., T. Stadler, S. Klopfstein, T. A. Heath, and F. Ronquist. 2016. 
#' Total-Evidence Dating under the Fossilized Birth-Death Process.
#' \emph{Systematic Biology} 65(2):228-249. 

#' @seealso 
#' \code{\link{createMrBayesConstraints}}, \code{\link{createMrBayesTipDatingNexus}}

#' @examples
#' 
#' # load retiolitid dataset
#' data(retiolitinae)
#' 
#' # uniform prior, with a 10 million year offset for
#' 	# the expected tree age from the earliest first appearance
#' 
#' createMrBayesTipCalibrations(tipTimes=retioRanges, whichAppearance="first",
#' 	ageCalibrationType="uniformRange", treeAgeOffset=10)
#' 
#' # fixed prior, at the earliest bound for the first appearance
#' 
#' createMrBayesTipCalibrations(tipTimes=retioRanges, whichAppearance="first",
#' 	ageCalibrationType="fixedDateEarlier", treeAgeOffset=10)
#' 
#' # fixed prior, sampled from between the bounds on the last appearance
#' 	# you should probably never do this, fyi
#' 
#' createMrBayesTipCalibrations(tipTimes=retioRanges, whichAppearance="first",
#' 	ageCalibrationType="fixedDateRandom", treeAgeOffset=10)
#' 
#' 
#' \dontrun{
#' 
#' createMrBayesTipCalibrations(tipTimes=retioRanges, whichAppearance="first",
#' 	ageCalibrationType="uniformRange", treeAgeOffset=10, file="tipCalibrations.txt")
#' 
#' }
#' 


#' @name createMrBayesTipCalibrations
#' @rdname createMrBayesTipCalibrations
#' @export
createMrBayesTipCalibrations<-function(tipTimes,
	ageCalibrationType,whichAppearance="first",
	treeAgeOffset,minTreeAge=NULL,file=NULL){
	#
	#
	if(length(ageCalibrationType)!=1){
		stop("argument ageCalibrationType must be of length 1")
		}
	if(length(whichAppearance)!=1){
		stop("argument whichAppearance must be of length 1")
		}
	if(all(whichAppearance!=c("first","last"))){
		stop("argument whichAppearance must be one of 'first' or 'last'")
		}
	if(all(ageCalibrationType!=c(
		"fixedDateEarlier","fixedDateLatter",
		"fixedDateRandom","uniformRange"))){
			stop('argument ageCalibrationType must be one of 
				"fixedDateEarlier", "fixedDateLatter", "fixedDateRandom", or "uniformRange"')
		}
	if(length(treeAgeOffset)!=1){
		stop("treeAgeOffset must be of length 1")
		}
	if(!is.null(minTreeAge)){
		if(length(minTreeAge)!=1){
			stop("minTreeAge must be of length 1")
			}
		}
	#
	############################################
	# format tipTimes
	#
	#if list of two (i.e. timeList), convert to four-date format
	if(is.list(tipTimes) & length(tipTimes)==2){
		tipTimes<-timeList2fourDate(tipTimes)
		}
	# checks for insanity in tipTimes
	if(!is.matrix(tipTimes)){
		stop("tipTimes must be of type matrix, or a list with a length of 2")
		}
	# must be a table with 1, 2 or 4 columns
	if(all(ncol(tipTimes)!=c(1,2,4))){
		stop("tipTimes must have 1 or 2 or 4 columns")
		}
	# Before we go any further, let's preserve the taxon names
	taxonNames<-rownames(tipTimes)
	#
	# check for things the user is unlikely to want to do
	if(ncol(tipTimes)==1){
		if(ageCalibrationType!="fixedDateEarlier"){
			stop("You appear to be supplying a single point occurrence per taxon.
				There isn't any uncertainty or upper bounds on ages, so 
				ageCalibrationType should be set to 'fixedDateEarlier'")
			}
		}
	####################################################################
	# fix so four columns
	# if two columns, then make four-date by repeating
	if(ncol(tipTimes)==2){
		tipTimes<-cbind(tipTimes,tipTimes)
		}
	#if a single column of point ages, then repeat four times
	if(ncol(tipTimes)==1){
		tipTimes<-cbind(tipTimes,tipTimes,tipTimes,tipTimes)
		}
	# choose either the first or last appearance times
		# (these will likely often be identical)
	if(whichAppearance=="first"){
		tipTimes<-tipTimes[,1:2]
		}
	if(whichAppearance=="last"){
		tipTimes<-tipTimes[,3:4]
		}
	#check
	if(ncol(tipTimes)!=2){
		stop("Weird data format for tipTimes")
		}
	#
	#########################################################
	# select times
	# for converting to fixed or uniform age ranges
	#
	if(ageCalibrationType=="fixedDateEarlier"){
		# use lower bound for selected age of appearance
		tipTimes<-tipTimes[,1,drop=FALSE]
		timeType<-"fixed"
		}
	if(ageCalibrationType=="fixedDateLatter"){	
		# use upper bound for selected age of appearance
		tipTimes<-tipTimes[,2,drop=FALSE]
		timeType<-"fixed"
		}
	if(ageCalibrationType=="fixedDateRandom"){	
		# random drawn from a uniform distribution
		tipTimes<-t(t(apply(tipTimes,1,function(x) runif(1,x[2],x[1]))))
		timeType<-"fixed"
		}
	if(ageCalibrationType=="uniformRange"){
		# if uniform, don't need to edit tipTimes
		timeType<-"uniform"
		}
	# check
	if(all(timeType!=c("fixed","uniform"))){
		stop("Problem when selecting ages. ageCalibrationType argument incorrect?")
		}
	# 
	##############################################
	# start writing MrBayes block
	if(timeType=="fixed"){
		#format fixed age script - single age per taxon
		dateBlock<-sapply(1:nrow(tipTimes),function(i)
			paste0("calibrate ",rownames(tipTimes)[i],
				" = fixed (",tipTimes[i],");"))
		}
	# use upper and lower bounds of selected age
		# of appearance to place uniform prior on tip age
	if(timeType=="uniform"){
		# format uniform age block - two ages per taxon
		dateBlock<-sapply(1:nrow(tipTimes),function(i)
			paste0("calibrate ",rownames(tipTimes)[i],
				" = uniform (",tipTimes[i,2],
				", ",tipTimes[i,1],");"))
		}
	#####################################################
	#need to create tree age prior
	# get minimum age of tips
	minTipAge<-max(tipTimes)
	if(is.null(minTreeAge)){
		minTreeAge<-minTipAge
	}else{
		# make sure minimum tree age is less than oldest tip
		if((minTreeAge*1.0001)<minTipAge){
			stop("User given minTreeAge is younger than the oldest tip age")
			}
		}
	#
	# use offset to calculate mean tree age
	meanTreeAge<-minTreeAge+treeAgeOffset
	#
	# write tree age prior command
	treeAgeBlock<-paste0("prset treeagepr = offsetexp(",
		minTreeAge,", ",meanTreeAge,");")
	########################################################
	# create final block for output
	#
	finalBlock<-c(dateBlock,"",treeAgeBlock)
	#
	if(!is.null(file)){
		write(finalBlock,file)
	}else{
		return(finalBlock)
		}
	}

