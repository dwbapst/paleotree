#' Calculating Diversity Curves Across Multiple Datasets
#' 
#' Calculates multiple diversity curves from a list of datasets of taxic ranges
#' and/or phylogenetic trees, for the same intervals, for all the individual
#' datasets. A median curve with 95 percent quantile bounds is also calculated
#' and plotted for each interval.
#' 
#' @details
#' This function is essentially a wrapper for the individual diversity curve
#' functions included in paleotree. \code{multiDiv} will intuitively decide whether
#' input datasets are continuous-time taxic ranges, discrete-time (binned
#' interval) taxic ranges or phylogenetic trees, as long as they are formatted
#' as required by the respective diversity curve functions. A list that
#' contains a mix of data types is entirely acceptable.
#' A list of matrices output from \code{fossilRecord2fossilTaxa},
#' via simulation with \code{simFossilRecord} is allowable,
#' and treated as input for \code{taxicDivCont}.
#' Data of an unknown type gives back an error.
#' 
#' The argument \code{split.int} splits intervals, if and \emph{only} if discrete interval
#' time data is included among the datasets. See the help file for \code{taxicDivDisc}
#' to see an explanation of why \code{split.int = TRUE} by default is probably a good
#' thing.
#' 
#' As with many functions in the \code{paleotree} library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' The 'averaged' curve is actually the median rather than the mean as
#' diversity counts are often highly skewed (in this author's experience).
#' 
#' The shaded certainty region around the median curve is the two-tailed 95
#' percent lower and upper quantiles, calculated from the observed data. It is
#' not a true probabilisitic confidence interval, as it has no relationship to
#' the standard error.

#' @aliases multiDiv plotMultiDiv

#' @inheritParams DiversityCurves

#' @param data A list where each element is a dataset, formatted to be input in
#' one of the diversity curve functions listed in \code{\link{DiversityCurves}}.

#' @param plot If \code{TRUE}, the median diversity curve is plotted.

#' @param results The output of a previous run of \code{multiDiv} for replotting.

#' @param plotMultCurves If \code{TRUE}, each individual diversity curve is plotted
#' rather than the median diversity curve and 95 percent quantiles. 
#' \code{plotMultCurves = FALSE} by default.

#' @param yAxisLims Limits for the y (i.e. richness) axis on the plotted diversity curves. 
#' Only affects plotting. Given as either \code{NULL} (the default) or as a vector of
#' length two as for \code{xlim} in the basic R function \code{plot}. 
#' Time axes will be plotted \emph{exactly} to these values. 
#' The minimum value must be more than 1 if \code{plotLogRich = TRUE}.

#' @param multRainbow If \code{TRUE} and \code{plotMultCurves = TRUE}, each line is
#' plotted as a different, randomized color using the function \code{rainbow}. If
#' \code{FALSE}, each line is plotted as a black line. This argument is ignored if
#' \code{divPalette} is supplied.

#' @param divPalette Can be used so users can pass a vector of chosen color
#' identifiers for each diversity curve in \code{data} which will take precedence
#' over \code{multRainbow}. Must be the same length as the number of diversity curves
#' supplied.

#' @param divLineType Used to determine line type (\code{lty}) of the
#' diversity curves plotted when \code{plotMultCurves  = } \code{TRUE}. Default
#' is \code{lty = 1} for all curves. Must be either length of 1 or 
#' exact length as number of diversity curves.

#' @param main The main label for the figure.

#' @return A list composed of three elements will be invisibly returned:
#' \item{int.times}{A two column matrix giving interval start and end times}
#' \item{div}{A matrix of measured diversities in particular intervals by rows,
#' with each column representing a different dataset included in the input}
#' \item{median.curve}{A three column matrix, where the first column is the
#' calculated median curve and the second and third columns are the 95 percent
#' quantile upper and lower bounds}

#' @seealso The diversity curve functions used include: \code{\link{phyloDiv}},
#' \code{\link{taxicDivCont}} and \code{\link{taxicDivDisc}}.
#' 
#' Also see the function \code{LTT.average.root} in the package \code{TreeSim}, which
#' calculates an average LTT curve for multiple phylogenies, the functions
#' \code{mltt.plot} in ape and \code{ltt} in \code{phytools}.

#' @examples
#' # let's look at this function
#'     # with some birth-death simulations
#' 
#' set.seed(444)
#' 
#' # multiDiv can take output from simFossilRecord
#'     # via fossilRecord2fossilTaxa
#' 
#' # what do many simulations run under some set of
#'     # conditions 'look' like on average?
#' set.seed(444)
#' records <- simFossilRecord(
#'     p = 0.1, 
#'     q = 0.1, 
#'     nruns = 10,
#'     totalTime = 30, 
#'     plot = TRUE
#'     )
#' 
#' taxa <- lapply(records, fossilRecord2fossilTaxa)
#' 
#' multiDiv(taxa)
#' # increasing cone of diversity! 
#' 
#' # Its even better on a log scale:
#' multiDiv(taxa, plotLogRich = TRUE)
#' 
#' #######################################
#' # pure-birth example with simFossilRecord
#' # note that conditioning is tricky
#' 
#' set.seed(444)
#' recordsPB <- simFossilRecord(
#'     p = 0.1, 
#'     q = 0, 
#'     nruns = 10,
#'     totalTime = 30,
#'     plot = TRUE
#'     )
#'     
#' taxaPB <- lapply(recordsPB, fossilRecord2fossilTaxa)
#' multiDiv(taxaPB, plotLogRich = TRUE)
#' 
#' #compare many discrete diversity curves
#' discreteRanges <- lapply(taxaPB, function(x)
#'     binTimeData(
#'         sampleRanges(x, 
#'             r = 0.5,
#'             min.taxa = 1
#'             ),
#'         int.length = 7)
#'     )
#' 
#' multiDiv(discreteRanges)
#' 
#' #########################################
#' # plotting a multi-diversity curve for
#'    # a sample of stochastic dated trees
#' 
#' record <- simFossilRecord(
#'     p = 0.1, q = 0.1, 
#'     nruns = 1,
#'     nTotalTaxa = c(30,40), 
#'     nExtant = 0)
#'     
#' taxa <- fossilRecord2fossilTaxa(record)
#' rangesCont <- sampleRanges(taxa, r = 0.5)
#' rangesDisc <- binTimeData(rangesCont,
#'     int.length = 1)
#' # get the cladogram
#' cladogram <- taxa2cladogram(taxa, plot = TRUE)
#' 
#' #using multiDiv with samples of trees
#' ttrees <- timePaleoPhy(
#'     cladogram, 
#'     rangesCont, 
#'     type = "basic",
#'     randres = TRUE, 
#'     ntrees = 10, 
#'     add.term = TRUE
#'     )
#'     
#' multiDiv(ttrees)
#' 
#' # uncertainty in diversity history is solely due to 
#'     # the random resolution of polytomies
#' 
#' ######################################################### 
#' 
#' #using multiDiv to compare very different data types:
#'    # continuous ranges, discrete ranges, dated tree
#' 
#' # get a single dated tree
#' ttree <- timePaleoPhy(
#'     cladogram, 
#'     rangesCont, 
#'     type = "basic", 
#'     add.term = TRUE, 
#'     plot = FALSE
#'     )
#'     
#' # put them altogether in a list
#' input <- list(rangesCont, rangesDisc, ttree)
#' 
#' multiDiv(input, plot = TRUE)
#' 
#' # what happens if we use fixed interval times?
#' multiDiv(input, 
#'     int.times = rangesDisc[[1]], 
#'     plot = TRUE)
#' 
#' layout(1)
#' 

#' @rdname multiDiv
#' @export
multiDiv <- function(
		data,
		int.length = 1,
		plot = TRUE,
		split.int = TRUE,
		drop.ZLB = TRUE,
		drop.cryptic = FALSE,
		extant.adjust = 0.01,
		plotLogRich = FALSE,
		yAxisLims = NULL,
		timelims = NULL,
		int.times = NULL,
		plotMultCurves = FALSE,
		multRainbow = TRUE,
		divPalette = NULL,
		divLineType = 1,
		main = NULL
		){
	###############################################
	#lines up a bunch of taxic or phylo objects and calculates diversity curves simulataneously
		#across all their objects; intuits the type of object without being told
		#it also calculates a "average" median curve and 95% quantile intervals
	#input is a list of dicrete interval or continuous time taxic data or a timetree
		#as in the respective functions
	#output is a list with third objects
		#the first object is a 2-column matrix with interval starts and ends
		#the second object is a matrix 
			#with the measured diversity for all the objects as columns, intervals as rows
	#3rd object consists of a 3col matrix of information related to median curve
		#first column is a per-interval median of the combined diversity curves
		#second and third columns are 95% quantile intervals on that median
	#int.length = 1;plot = TRUE;split.int = TRUE;drop.ZLB = TRUE;drop.cryptic = FALSE;plotLogRich = FALSE;timeLims = NULL
	#plotMultCurves = FALSE;multRainbow = TRUE;divPalette = NULL
	#require(ape)
	#
  # identify class of data
  # data classes
  dataClassType <- sapply(data, function(x){
    classType <- NA
    if(inherits(x = x, what = "matrix")){
      classType <- 1
      }
    if(inherits(x = x, what = "list")){
      classType <- 2
      }  
    if(inherits(x = x, what = "phylo")){
      classType <- 3
      }
    if(is.na(classType)){
      stop("Data of unknown type - not a matrix, list or phylogeny?")
      }
    return(classType)
    })
  # 
	if(is.null(int.times)){
		tblen <- int.length
		#get max and min times for each type
		if(any(dataClassType == 1)){
			lims1 <- sapply(data[dataClassType == 1],function(x) 
				if(ncol(x) == 6){
					c(min(x[,3:4],na.rm = T),max(x[,3:4],na.rm = T))	
				}else{
					c(min(x,na.rm = TRUE),max(x,na.rm = TRUE))}
				)
		}else{lims1 <- NA}
		if(any(dataClassType == 2)){	
			for(i in which(dataClassType == 2)){
				data[[i]][[1]][data[[i]][[1]][,1] == 0,1] <- extant.adjust
				}
			lims2 <- sapply(data[dataClassType == 2],function(x) 
				c(min(x[[1]][max(x[[2]]),]),max(x[[1]][min(x[[2]]),])))
		}else{lims2 <- NA}
		if(any(dataClassType == 3)){
			lims3 <- numeric()
			for(i in which(dataClassType == 3)){
				ttree <- data[[i]]
				if(is.null(ttree$root.time)){
					ntime <- node.depth.edgelength(ttree)
					ntime <- max(ntime)-ntime
				}else{
					ntime <- node.depth.edgelength(ttree)
					ntime <- ttree$root.time-ntime
					}
				lims3 <- c(lims3,c(min(ntime),max(ntime)))
				}
		}else{lims3 <- NA}
		end <- min(c(lims1,lims2,lims3),na.rm = TRUE)
		start <- max(c(lims1,lims2,lims3),na.rm = TRUE)
		midtimes <- seq(start+2*tblen,end-2*tblen,by = -tblen)
		midtimes <- midtimes[midtimes>0]
		int.start <- midtimes+(tblen/2);int.end <- midtimes-(tblen/2)
		int.times <- cbind(int.start,int.end)
		if(split.int & any(dataClassType == 2)){
			# for every single discrete time dataset
		      # bins must be split at their boundaries
			splinters <- sapply(data[dclasss = 2],function(x) x[[1]][,1])
			mustSplit <- apply(int.times, 1,
				function(x) any(sapply(splinters,function(y) x[1]>y & x[2]<y))
				)
			if(any(mustSplit)){
				for(i in which(mustSplit)){
						splitter <- splinters[sapply(splinters[,1],
							function(y) int.times[i,1]>y & int.times[i,2]<y
							),1]
						#if(length(splitter)>1){stop("Splitter returning more than one value?!")}
						splitter <- c(int.times[i,1],splitter,int.times[i,2])
						int.times <- rbind(
							int.times,
							cbind(
								splitter[1:(length(splitter)-1)],
								splitter[2:length(splitter)]
								)
							)
					}
				int.times <- int.times[-which(mustSplit),]
				int.times <- int.times[order(-int.times[,1]),]
				}
			midtimes <- (int.start+int.end)/2		
			}
		}			
	div <- matrix(,nrow(int.times),1)
	for(i in 1:length(data)){
		if((dataClassType[i] == 1)){
			divs1 <- taxicDivCont(
				timeData = data[[i]],
				int.times = int.times,
				plot = FALSE,
				drop.cryptic = drop.cryptic
				)
			divs1 <- divs1[,3]
			div <- cbind(div,divs1)
			}
		if((dataClassType[i] == 2)){
			divs2 <- taxicDivDisc(
				timeList = data[[i]],
				int.times = int.times,
				split.int = FALSE,
				plot = FALSE
				)
			divs2 <- divs2[,3]
			div <- cbind(div,divs2)
			}
		if((dataClassType[i] == 3)){
			divs3 <- phyloDiv(
				tree = data[[i]],
				int.times = int.times,
				plot = FALSE,
				drop.ZLB = drop.ZLB
				)
			divs3 <- divs3[,3]
			div <- cbind(div,divs3)
			}
		}
	div <- div[,-1]
	colnames(div) <- paste("dataset",1:ncol(div),sep = "") 
	#get median curve
	median <- apply(div,1,median)
	#the low quantile
	q1 <- apply(div,1,quantile,probs = 0.025)	
	#the high quantile
	q2 <- apply(div,1,quantile,probs = 0.975)	
	median.curve <- cbind(
		median = median,
		low.95.quantile = q1,
		high.95.quantile = q2
		)
	#
	res <- list(
		int.times = int.times,
		div.counts = div,
		median.curve = median.curve
		)
	#
	if(plot){
		plotMultiDiv(
			res,
			plotLogRich = plotLogRich,
			timelims = timelims, 
			yAxisLims = yAxisLims,
			plotMultCurves = plotMultCurves,
			multRainbow = multRainbow,
			divPalette = divPalette,
			divLineType = divLineType,
			main = main
		)}
	#
	return(invisible(res))
	}

#' @rdname multiDiv
#' @export	
plotMultiDiv <- function(
      results, 
      plotLogRich = FALSE, 
      timelims = NULL, 
      yAxisLims = NULL, 
      plotMultCurves = FALSE, 
      multRainbow = TRUE,
      divPalette = NULL,
      divLineType = 1,
      main = NULL
      ){
  ####################
	if(!is.null(timelims)){
		if(length(timelims) != 2){
			stop("if given, timelims should be a vector of length 2")
			}
		}
	if(!is.null(yAxisLims)){
		if(length(yAxisLims) != 2){
			stop("if given, yAxisLims should be a vector of length 2")
			}
		if(plotLogRich & min(yAxisLims)<1){
			stop("if richness is plotted on a log scale, minimum axis value must be  >= 1")
			}
		}
	# plots the median diversity curve for a multiDiv() result
	int.start <- results[[1]][,1]
	int.end <- results[[1]][,2]
	times1 <- c(int.start,(int.end+((int.start-int.end)/10000)))
	if(plotMultCurves){
		if(is.null(main)){
			main <- "Multiple Diversity Curves"
		  }
	  # here's my diversity information
		divs <- results[[2]]	
		divs1 <- rbind(divs,divs)[order(times1),]
		times1 <- sort(times1)
		# set up the general plotting window
		if(plotLogRich){
			if(is.null(yAxisLims)){
				y_lim <- c(min(divs1[divs1 >= 1]),max(divs1[divs1 >= 1]))
			}else{
				y_lim <- yAxisLims
				}
			plot(times1[divs1[,1]>0],divs1[divs1[,1]>0,1],type = "n",
				ylim = y_lim,log = "y",
				xlim = if(is.null(timelims)){
						c(max(times1),max(0,min(times1)))
					}else{
						timelims
						},
				xaxs = if(is.null(timelims)){"r"}else{"i"},
				xlab = "Time (Before Present)",
				ylab = "Log Lineage/Taxic Richness",
				main = )
		}else{
			if(is.null(yAxisLims)){
				y_lim <- c(min(divs1),max(divs1))
			}else{
				y_lim <- yAxisLims
				}
			plot(times1,divs1[,1],type = "n",
				ylim = y_lim,
				xlim = if(is.null(timelims)){
						c(max(times1),max(0,min(times1)))
					}else{
						timelims
						},
				xaxs = if(is.null(timelims)){"r"}else{"i"},
				xlab = "Time (Before Present)",
				ylab = "Lineage/Taxic Richness",
				main = main)
			}
		if(is.null(divPalette)){
			if(multRainbow){
				divPalette <- sample(rainbow(ncol(divs1)))
			}else{
				divPalette <- rep(1,ncol(divs1))
				}
			}
		if(length(divLineType) != ncol(divs1)){
			divLineType <- rep(divLineType,ncol(divs1))
			}
		for(i in 1:ncol(divs1)){	#plot each line
			lines(times1,
				divs1[,i],
				lwd = 3,
				col = divPalette[i],
				lty = divLineType[i]
				)
			}
	}else{
		if(is.null(main)){
			main <- "Median Diversity Curve"
			}
		mdiv <- results[[3]]
		mdiv1 <- rbind(mdiv,mdiv)[order(times1),]
		times1 <- sort(times1)
		if(plotLogRich){
			mdiv1[mdiv1[,2]<1,2] <- 1
			mdiv1[mdiv1[,3]<1,3] <- 1
			if(is.null(yAxisLims)){
				y_lim <- c(min(mdiv1[mdiv1 >= 1]), max(mdiv1[mdiv1 >= 1]))
			}else{
				y_lim <- yAxisLims
				}
			plot(times1[mdiv1[,3]>0],mdiv1[mdiv1[,3]>0,3],type = "n",
				ylim = y_lim,log = "y",
				xlim = if(is.null(timelims)){
						c(max(times1),max(0,min(times1)))
					}else{
						timelims
						},
				xaxs = if(is.null(timelims)){"r"}else{"i"},
				xlab = "Time (Before Present)",
				ylab = "Log Lineage/Taxic Richness",
				main = main)
		}else{
			if(is.null(yAxisLims)){
				y_lim <- c(min(mdiv1),max(mdiv1))
			}else{
				y_lim <- yAxisLims
				}
			plot(times1,mdiv1[,3],type = "n",
				ylim = y_lim,
				xlim = if(is.null(timelims)){
						c(max(times1),max(0,min(times1)))
					}else{
						timelims
						},
				xaxs = if(is.null(timelims)){"r"}else{"i"},
				xlab = "Time (Before Present)",
				ylab = "Lineage/Taxic Richness",
				main = main)
			}
		polygon(
			c(times1,rev(times1)),
			c(mdiv1[,3],rev(mdiv1[,2])),
			col = "gray",
			border = NA
			)
		#
		lines(times1,mdiv1[,1],lwd = 3)
		#
		}
	}