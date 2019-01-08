





#' @par phylopicIDsPBDB ID numbers for images from Phylopic
#', as given by the Paleobiology Database's API under the output
#' \code{image_no} (given when \code{show = img}).

# note size is vertical, proportional to the space between tips
	# max horizontal size stops flat / long phylopics from becoming overly huge

# set x.lim so plot x limits is * (1 + extraMargin)
# where 1 is the tree height (effectively)

# noiseThreshold threshold for noise in the PNG from Phylopic to be
# treated as meaningless noise (i.e. a color that is effectively whitespace)
# and thus can be trimmed as margin to be trimmed by the function



plotPhylopicTreePBDB <- function(tree, 
		phylopicIDsPBDB = NULL, 
		size = 0.9,
		noiseThreshold = 0.1,
		extraMargin = 0.2,
		rescalePNG = TRUE,
		trimPNG = TRUE
		){		
	#########################################
	# uses calls to the Paleobiology Database's API
		# to construct a phylogeny with phylopics
			# (from Phylopic, duh)
		# as pictorial replacements for the tip labels
	############
	#	
	# require(png);require(RCurl);require(ape)
	source("D:\\dave\\workspace\\paleotree\\R\\plotPhylopicTreePBDB.R")
	#
	# check or obtain the phylopic IDs from PBDB
	phylopicIDsPBDB <- getPhyloPicIDNum(
		phylopicIDsPBDB = phylopicIDsPBDB,
		tree = tree)
	##############################################
	# plot a tree but with blank tip labels
		# set x.lim so plot x limits is * (1 + extraMargin)
		# where 1 is the tree height (effectively)
	# calculate new x.lim by
		# *not* plotting a tree
	outPlot <- plot.phylo(tree,
		plot=FALSE,
		show.tip.label=FALSE)
	old_xlim <- outPlot$x.lim[2]
	new_xlim <- old_xlim * (1 + extraMargin)
	par(new = TRUE)
	plot.phylo(tree,
		x.lim = new_xlim,
		show.tip.label = FALSE)
	##########################################
	# now get the last plotting environment
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	# get the device's aspect ratio
	devAspRatio <- dev.size()[1] / dev.size()[2]
	# get the plot's own aspect ratio
	plotAspRatio <- diff(lastPP$x.lim) / diff(lastPP$y.lim)
	# true aspect ratio is their product apparently
	plotAspRatio <- plotAspRatio * devAspRatio 
	#
	# pause 3 seconds so we don't spam the API
	Sys.sleep(3)
	for (i in 1:lastPP$Ntip){
		# GET IMAGE
		picPNG <- getPhyloPicPNG(picID = phylopicIDsPBDB[i], 
			noiseThreshold = noiseThreshold,
			rescalePNG = rescalePNG,
			trimPNG = trimPNG)
		#
		#########################################
		#get aspect ratio
			# ratio of # of pixel dimensions
		picAspRatio <- dim(picPNG)[1]/dim(picPNG)[2]
		#############################################
		# adjustment of size
			# need to modify size relative to aspect ratio
		if(picAspRatio < 1){
			# its flatish so correct it by aspect ratio
			picSize <- size * (picAspRatio^0.7)
		}else{
			# then its skinny, not flat
			# don't do anything
			picSize <- size 
			}
		###########################################
		# GET THE COORDINATES
		#
		# offset is size/2 by default
		offset <- size*0.9/plotAspRatio
		#
		#points(lastPP$xx,lastPP$yy)	
		x<-lastPP$xx[i]+offset
		y<-lastPP$yy[i]
		#
		##################################################
		#adjust the position of the sides for the image
		#
		xAdj <- (picSize /2) / (picAspRatio*plotAspRatio) 
		yAdj <- picSize /2
		#
		##################################################
		# plot the picPNG using graphics::rasterImage
		graphics::rasterImage(picPNG,
			xleft = x - xAdj ,
			ybottom = y - yAdj ,
			xright = x + xAdj ,
			ytop = y + yAdj,
			interpolate = TRUE
			)
		# cool
		}
	# what to return? nothing I guess
	}


getPhyloPicPNG<-function(picID, noiseThreshold=0.1,
		rescalePNG = TRUE, trimPNG = TRUE,
		plotComparison = FALSE){
	############################################
	# GET IMAGE
	# get the URL address for the pic via API
	apiPicURL <- paste0(
		"http://paleobiodb.org/data1.2/taxa/thumb.png?id=",
		picID)
	# pause 1 second
	Sys.sleep(runif(1,1,2))
	# get picPNG		
	picPNG <-  png::readPNG(RCurl::getURLContent(apiPicURL))
	############
	# phylopic PNG is on the fourth slice
	#  image(picPNG [,,4])
	########################
	# RESCALE PALETTE
	sliceOriginal <- picPNG[,,4]
	if(rescalePNG){
		#rescale pic so that min is 0 and max is 1
		picPNG[,,4] <- picPNG[,,4]-abs(min(picPNG[,,4]))
		picPNG <- picPNG/max(picPNG)
		}
	#
	###################
	# TRIM THE PHYLOPIC
		# lots of phylopics have contiguous whitespace at the top/bottom
	sliceContrasted <- picPNG[,,4]
	if(trimPNG){
		sliceContrasted[sliceContrasted  < noiseThreshold] <- 0
		# find all rows of the PNG from the top AND bottom
			# that are non-contiguous whitespace
			# these are to be SAVED
			# use 'contrasted' version
		saveRows <- ((cumsum(apply(sliceContrasted , 1, sum)) > 0)
			 & rev(cumsum(rev(apply(sliceContrasted , 1, sum))) > 0))
		#
		# turns out lots phylopics also have whitespace on their left/right
			# need to trim this too
		# find all COLUMNS of the PNG from the RIGHT AND LEFT
			# that are non-contiguous whitespace
			# these are to be SAVED
		saveCols <- ((cumsum(apply(sliceContrasted , 2, sum)) > 0)
			 & rev(cumsum(rev(apply(sliceContrasted , 2, sum))) > 0))
		#
		# remove whitespace from all slices of the array
		picPNG <- picPNG[saveRows,saveCols,]
		}
	##############
	#
	if(plotComparison){
		# plots a comparison of three images
			# as a diagnostic mode
		layout(1:3)
		par(mar=c(0,0,0,0))
		image(sliceOriginal)
		image(sliceContrasted)
		image(picPNG [,,4])
		}
	return(picPNG)
	}



getPhyloPicIDNum <- function(phylopicIDsPBDB, tree){
	# check or obtain the phylopic IDs from PBDB
	#
	if(is.null(phylopicIDsPBDB)){
		# get image ID numbers using PBDB API calls for each
			# tip taxon in the tree using the tip labels
		tiptaxa <- paste0(tree$tip.label,collapse = ",")
		apiAddressTaxa <- paste0(
			"http://paleobiodb.org/data1.2/taxa/list.txt?name=",taxa,
			"&rel=exact&show=img"
			)	
		# call PBDB API
		tiptaxaData <- read.csv(apiAddressTaxa,
			stringsAsFactors = FALSE)
		# get the image IDs
		phylopicIDsPBDB<- tiptaxaData $image_no[
			match(tree$tip.label, tiptaxaData $taxon_name)]
		names(phylopicIDsPBDB) <- tree$tip.label
	}else{
		# CHECKS
		# does the provided number of IDs equal the number
			# of tips?
		if(length(phylopicIDsPBDB) != Ntip(tree)){
			#print(tree)
			#print(phylopicIDsPBDB)
			stop(paste0("If provided, phylopicIDsPBDB must be have the same",
				" length as the number of tips on the input tree"))
			}
		# make sure its a character vector
		if(!is.character(phylopicIDsPBDB)){
			phylopicIDsPBDB <- as.character(phylopicIDsPBDB)
			if(!is.character(phylopicIDsPBDB)){		
				stop("Cannot coerce phylopicIDsPBDB to be type character")
				}
			}
		# Are there NAs or blanks?
		if(any(is.na(phylopicIDsPBDB))){
			stop("If provided, phylopicIDsPBDB cannot contain any NA values")
			}
		if(any(phylopicIDsPBDB == "")){
			stop("If provided, phylopicIDsPBDB cannot contain any blank values")
			}
		# make sure user provided IDs are labeled
		if(length(names(phylopicIDsPBDB)) != Ntip(tree)){
			print(length(names(phylopicIDsPBDB)))
			print(phylopicIDsPBDB)
			stop(paste0("If provided, phylopicIDsPBDB must have",
				" labels matching the tree's tip labels"))	
			}
		}
	# make sure phylopicIDsPBDB is sorted relative to tip order
	phylopicIDsPBDB <- phylopicIDsPBDB[
		match(tree$tip.label, names(phylopicIDsPBDB))]
	return(phylopicIDsPBDB)
	}
