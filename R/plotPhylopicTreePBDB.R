#' Plot a Phylogeny with Organismal Silhouettes from PhyloPic, Called Via the Paleobiology Database API
#'
#' This function will take a phylogeny, preferably a taxonomy-tree
#' created from classification information and/or parent-child taxon
#' information pulled from the Paleobiology Database via function
#' \code{\link{makePBDBtaxonTree}}, and use the
#' Paleobiology Database's API to plot silhouettes of each given tip taxon
#' in replacement of their normal tip labels.


#' @details
#' This function preferably will pull the identifiers for which images are to be associated with the tip taxa from \code{taxaDataPBDB$image_no}. By default, \code{taxaDataPBDB} itself is assumed to be an element of \code{tree} named \code{tree$taxaData}, as the PBDB data table used to construct the tree is appended to the output tree when \code{\link{makePBDBtaxonTree}} is used to construct a taxonomy-tree. If the functions listed in \code{\link{getTaxaDataPBDB}} are used to obtain the taxonomic data, this table will include the \code{image_no} variable, which is the image identifier numbers needed to call PNGs from the Paleobiology Database API. If \code{taxaDataPBDB} isn't provided, either by the user directly, or as an element of \code{tree}.




#' @param tree A phylogeny of class \code{phylo} which will be
#' plotted, with the terminal tip taxa replaced by silhouettes.
#' The tree will be plotted with edge lengths.

#' @param taxaDataPBDB  A \code{data.frame} of taxonomic data from
#' the Paleobiology Database containing an \code{$image_no} variable,
#' as returned when \code{show = "img"} is used. See \emph{Details}.

#' @param size The default is \code{size = 0.9}.

#' @param noiseThreshold A threshold for noise in the PNG from Phylopic
#' to be treated as meaningless noise (i.e. a color that is effectively
#' whitespace) and thus can be trimmed as empty margin which can be
#' trimmed before the silhouette is plotted. The units for this argument
#' are on a scale from 0 to 1, with 0 being true white space, and values
#' between 0 and 0.5 representing colors closer to whitespace than true
#' black. The default is \code{noiseThreshold = 0.1}.



#' @param extraMargin The default is \code{extraMargin = 0.2}.

#' @param rescalePNG If \code{TRUE} (the default),

#' @param trimPNG If \code{TRUE} (the default),

#' @param makeMonochrome If \code{TRUE}, PhyloPic silhouettes are
#' forced to be purely monochrome black-and-white, with no gray
#' scale. Most of the silhouettes are binary black-and-white already
#' but some aren't, but those gray-scale values (sometimes?) seem
#' to exist to indicate very fine features. However, maybe an image
#' is far too much gray-scale, in which case users can try this
#' option to force all silhouettes to be monochrome.
#' The default is \code{FALSE}.

#' @param cacheDir If not \code{NULL}, first look here for a cached
#' version of the images. This makes loading faster.
#' The default is \code{NULL}.


# @param phylopicIDsPBDB ID numbers for images from Phylopic,
# as given by the Paleobiology Database's API under the output
# \code{image_no} (given when \code{show = img}).

# note size is vertical, proportional to the space between tips
	# max horizontal size stops flat / long phylopics from becoming overly huge

# set x.lim so plot x limits is * (1 + extraMargin)
# where 1 is the tree height (effectively)

# noiseThreshold threshold for noise in the PNG from Phylopic to be
# treated as meaningless noise (i.e. a color that is effectively whitespace)
# and thus can be trimmed as margin to be trimmed by the function


#' @return
#' This function returns nothing - it just plots the input tree with
#' silhouettes from the PhyloPic database.

#' @seealso
#' See \code{\link{getTaxaDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhylopicTreePBDB}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.
#'

#' @examples
#'
#' \donttest{
#'
#'
#' # now plot with phylopic images
#' plotPhylopicTreePBDB(tree = tree,
#' 	taxaDataPBDB = taxaData)
#'
#' }
#'


#' @name plotPhylopicTreePBDB
#' @rdname plotPhylopicTreePBDB
#' @export
plotPhylopicTreePBDB <- function(
		tree,
		taxaDataPBDB = tree$taxaDataPBDB,
		# phylopicIDsPBDB = NULL,
		size = 0.9,
		noiseThreshold = 0.1,
		extraMargin = 0.2,
		rescalePNG = TRUE,
		trimPNG = TRUE,
		makeMonochrome = FALSE,
		cacheDir = NULL
		){
	#########################################
	# uses calls to the Paleobiology Database's API
		# to construct a phylogeny with phylopics
			# (from Phylopic, duh)
		# as pictorial replacements for the tip labels
	###############################################
	# check or obtain the phylopic IDs from PBDB
	phylopicIDsPBDB <- getPhyloPicIDNum(
		taxaData = taxaDataPBDB,
		#phylopicIDsPBDB = phylopicIDsPBDB,
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
	devAspRatio <- grDevices::dev.size()[1] / grDevices::dev.size()[2]
	# get the plot's own aspect ratio
	plotAspRatio <- diff(lastPP$x.lim) / diff(lastPP$y.lim)
	# true aspect ratio is their product apparently
	plotAspRatio <- plotAspRatio / devAspRatio
	#
	# pause 3 seconds so we don't spam the API
	if(is.null(cacheDir)) {
		Sys.sleep(3)
	}
	for (i in 1:lastPP$Ntip){
		# GET IMAGE
		picPNG <- getPhyloPicPNG(picID = phylopicIDsPBDB[i],
			noiseThreshold = noiseThreshold,
			rescalePNG = rescalePNG,
			makeMonochrome = makeMonochrome,
			trimPNG = trimPNG,
		  cacheDir = cacheDir)
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
		offset <- size*0.9* plotAspRatio
		#
		#points(lastPP$xx,lastPP$yy)
		x<-lastPP$xx[i]+offset
		y<-lastPP$yy[i]
		#
		##################################################
		#adjust the position of the sides for the image
		#
		xAdj <- (picSize/2) * (plotAspRatio/picAspRatio)
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

#' Cache phylopic images for a tree
#'
#' This function takes a tree and saves phylopic thumbnails locally,
#' in order to dramatically speed up later rendering.
#'

#' @param tree A phylogeny of class \code{phylo} which will be
#' plotted, with the terminal tip taxa replaced by silhouettes.
#' The tree will be plotted with edge lengths.

#' @param taxaDataPBDB  A \code{data.frame} of taxonomic data from
#' the Paleobiology Database containing an \code{$image_no} variable,
#' as returned when \code{show = "img"} is used. See \emph{Details}.

#' @param cacheDir Where to save the output PNGs

#' @export

cachePhyloPicPNG <- function(
	tree, taxaDataPBDB= tree$taxaDataPBDB, cacheDir
){
	ids <- getPhyloPicIDNum(taxaData=taxaDataPBDB, tree=tree)
	for(i in seq_along(ids)) {
		picPNG <- getPhyloPicPNG(ids[i])
		magick::image_write(picPNG, path=file.path(cacheDir, paste0(ids[i], ".png")), format="png")
	}
}

getPhyloPicPNG<-function(
		picID, noiseThreshold = 0.1,
		rescalePNG = TRUE, trimPNG = TRUE,
		makeMonochrome = FALSE,
		plotComparison = FALSE,
	  cacheDir = NULL){
	############################################
	#
	# require(png);require(RCurl)
	# png::readPNG RCurl::getURLContent
	#
	# GET IMAGE
	#
	# First try to get a cached version
	picPNG <- NULL
	if(!is.null(cacheDir)) {
		try(picPNG <- png::readPNG(file.path(cacheDir, paste0(picID, ".png"))))
	}
	if(is.null(picPNG)) { # we have failed to load cached version
		# get the URL address for the pic via API
		apiPicURL <- paste0(
			"http://paleobiodb.org/data1.2/taxa/thumb.png?id=",
			picID)
		# pause 1 second
		Sys.sleep(runif(1,1,2))
		# get picPNG
		picPNG <-  png::readPNG(RCurl::getURLContent(apiPicURL))
	}
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
		#picPNG <- picPNG^0.75
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
	if(makeMonochrome){
		picPNG <- picPNG^0.001
		}
	##############
	#
	if(plotComparison){
		# plots a comparison of three images
			# as a diagnostic mode
		layout(1:3)
		par(mar=c(0,0,0,0))
		graphics::image(sliceOriginal)
		graphics::image(sliceContrasted)
		graphics::image(picPNG [,,4])
		}
	return(picPNG)
	}



getPhyloPicIDNum <- function(taxaData, tree){
	# check or obtain the phylopic IDs from PBDB
	#
	if(is.null(taxaData)){
		# get image ID numbers using PBDB API calls for each
			# tip taxon in the tree using the tip labels
		tiptaxa <- paste0(tree$tip.label,
			collapse = ",")
		apiAddressTaxa <- paste0(
			"http://paleobiodb.org/data1.2/taxa/list.txt?name=",
			tiptaxa, "&rel=exact&show=img"
			)
		# call PBDB API
		tiptaxaData <- read.csv(apiAddressTaxa,
			stringsAsFactors = FALSE)
		# get the image IDs
		phylopicIDsPBDB<- tiptaxaData $image_no[
			match(tree$tip.label, tiptaxaData $taxon_name)]
		names(phylopicIDsPBDB) <- tree$tip.label
	}else{
		#
		# get phylo pic IDs and label with tip labels
		phylopicIDsPBDB <- as.character(taxaData$image_no[
			match(tree$tip.label, taxaData$taxon_name)])
		names(phylopicIDsPBDB) <- tree$tip.label
		#
		# CHECKS
		# does the provided number of IDs equal the number
			# of tips?
		if(length(phylopicIDsPBDB) != Ntip(tree)){
			#print(tree)
			#print(phylopicIDsPBDB)
			stop(paste0(
				"If provided, phylopicIDsPBDB must be have the same",
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
		match(tree$tip.label, names(phylopicIDsPBDB))
		]
	return(phylopicIDsPBDB)
	}
