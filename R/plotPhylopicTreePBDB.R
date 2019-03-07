#' Plot a Phylogeny with Organismal Silhouettes from PhyloPic, Called Via the Paleobiology Database API
#' 
#' This function will take a phylogeny, preferably a taxonomy-tree
#' created from classification information and/or parent-child taxon
#' information pulled from the Paleobiology Database via function
#' \code{\link{makePBDBtaxonTree}}, and use the
#' Paleobiology Database's API to plot silhouettes of each given tip taxon
#' in replacement of their normal tip labels.

#' @details
#' This function preferably will pull the identifiers for which images are to
#' be associated with the tip taxa from \code{taxaDataPBDB$image_no}. By default,
#' \code{taxaDataPBDB} itself is assumed to be an element of \code{tree} named
#' \code{tree$taxaData}, as the PBDB data table used to construct the tree is
#' appended to the output tree when \code{\link{makePBDBtaxonTree}} is used to
#' construct a taxonomy-tree. If the functions listed in \code{\link{getTaxaDataPBDB}}
#' are used to obtain the taxonomic data, this table will include the \code{image_no}
#' variable, which is the image identifier numbers needed to call PNGs from the
#' Paleobiology Database API. If \code{taxaDataPBDB} isn't provided, either by
#' the user directly, or as an element of \code{tree}. 

#' @param tree A phylogeny of class \code{phylo} which will be
#' plotted, with the terminal tip taxa replaced by silhouettes.
#' The tree will be plotted with edge lengths.

#' @param taxaDataPBDB  A \code{data.frame} of taxonomic data from
#' the Paleobiology Database containing an \code{$image_no} variable,
#' as returned when \code{show = "img"} is used. See \emph{Details}.

#' @param sizeScale The default is \code{sizeScale = 0.9}.

#' @param noiseThreshold A threshold for noise in the PNG from Phylopic
#' to be treated as meaningless noise (i.e. a color that is effectively
#' whitespace) and thus can be trimmed as empty margin which can be
#' trimmed before the silhouette is plotted. The units for this argument
#' are on a scale from 0 to 1, with 0 being true white space, and values
#' between 0 and 0.5 representing colors closer to whitespace than true
#' black. The default is \code{noiseThreshold = 0.1}.

# noiseThreshold threshold for noise in the PNG from Phylopic to be
# treated as meaningless noise (i.e. a color that is effectively whitespace)
# and thus can be trimmed as margin to be trimmed by the function

#' @param removeSurroundingMargin This argument controls the \code{no.margin} argument
#' in the function \code{plot.phylo}, which controls whether a (very large) margin is 
#' placed around the plotted tree, or not. By default, \code{plotPhylopicTreePBDB} will
#' suppress that margin, so that the plotted tree goes (very nearly) to the edges
#' of the plotting area.

#' @param extraMargin How much extra margin should be added to the side of the graph
#' to which PhyloPics are being added? This value is 0.08 by default, which works well
#' if the margin surrounding the entire plot is suppressed via argument
#' \code{removeSurroundingMargin = TRUE}. If the surrounding margin is not suppressed,
#' plots look mostly okay if users change to \code{extraMargin = 2}. Obviously this value
#' should be tweaked for every tree, size of plot and aspect ratio for maximum clarity. 

#' @param rescalePNG If \code{TRUE} (the default), the downloaded PhyloPic 
#' has its color values rebalanced to go from the most extreme white
#' to the most extreme black. Some (especially PBDB's versions) have varying
#' levels of gray compression-related artifacts and may not be properly
#' on a black-to-white scale.

#' @param trimPNG If \code{TRUE} (the default), the PhyloPic PNG is trimmed 
#' to remove extraneous whitespace from the top and bottom, before rescaling 
#' of the color values of the PNG.

#' @param makeMonochrome If \code{TRUE}, PhyloPic silhouettes are
#' forced to be purely monochrome black-and-white, with no gray
#' scale. Most of the silhouettes are binary black-and-white already
#' but some aren't, but those gray-scale values (sometimes?) seem
#' to exist to indicate very fine features. However, maybe an image
#' is far too much gray-scale, in which case users can try this
#' option to force all silhouettes to be monochrome.
#' The default is \code{FALSE}.

# @param phylopicIDsPBDB ID numbers for images from Phylopic,
# as given by the Paleobiology Database's API under the output
# \code{image_no} (given when \code{show = img}).

#' @param cacheDir If not \code{NULL}, this value is used as 
#' the name of a sub-directory of the working directory for which to look for
#' (or store) cached versions of PhyloPic PNGs to save on processing speed
#' and the need to pull the images from an external PNG.
#' If \code{NULL}, then cached images will not be checked for, and images downloaded will not be cached.
#' The default is 

#' @param cacheImage If \code{TRUE} (the default), images downloaded from the Paleobiology
#' Database and/or the PhyloPic Database will be cached to save on processing speed
#' and avoid the need to pull the images from an external PNG.

#' @param taxaColor Controls the color of plotted PhyloPics. Can either be \code{NULL}
#' (the default, all taxa will be plotted as black), or a character vector that is either
#' length 1, or the same length as the number of taxa. If \code{taxaColor} is length 1,
#' then the value is either interpreted as matching a tip label (in which case, the
#' named taxon will be highlighted in bright red), or as a color, which all PhyloPics
#' will then be plotted as that color. If the vector is the same length as the number
#' of taxa on \code{tree}, each value should be a character value of a
#' named color in base R, allowing user control over each PhyloPic
#' individually.

#' @param transparency A numeric value between 0 and 1, either length 1, or the same
#' length as the number of tips on \code{tree}. This indicates the transparency of
#' either all the plotted PhyloPics, or allows user control over each PhyloPic
#' individually. The default is 1, which represents maximum opaqueness,
#' applied to all PhyloPics.
	
#' @param ... Additional arguments, passed to
#' \code{plot.phylo} for plotting of the tree. These
#' additional arguments may be passed to \code{plot},
#' and from there to \code{plot}.

#' @return
#' Thisfunction silently returns the positions for elements in the tree

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

# note sizeScale is vertical, proportional to the space between tips
	# max horizontal sizeScale stops flat / long PhyloPics from becoming overly huge

# set x.lim so plot x limits is * (1 + extraMargin)
# where 1 is the tree height (effectively)



#' @name plotPhylopicTreePBDB
#' @rdname plotPhylopicTreePBDB
#' @export
plotPhylopicTreePBDB <- function(
		tree, 
		taxaDataPBDB = tree$taxaDataPBDB,
		# phylopicIDsPBDB = NULL, 
		#######################
		sizeScale = 0.9,
		removeSurroundingMargin = TRUE,
		extraMargin = 0.08,
		###########################
		taxaColor = NULL,
		transparency = 1,
		######################
		cacheDir = "cachedPhyloPicPNGs",
		cacheImage = TRUE,		
		##########################
		noiseThreshold = 0.1,
		rescalePNG = TRUE,
		trimPNG = TRUE,
		makeMonochrome = FALSE,
		...
		){		
	#########################################
	# uses calls to the Paleobiology Database's API
		# to construct a phylogeny with PhyloPics
			# (from Phylopic, duh)
		# as pictorial replacements for the tip labels
	###############################################
	# check or obtain the phylopic IDs from PBDB
	phylopicIDsPBDB <- getPhyloPicIDNumFromPBDB(
		taxaData = taxaDataPBDB,
		tree = tree)
	###############################################
	# determine colors for every taxon using taxaColor
	taxaColor <- matchTaxaColor(
		taxaColorOld = taxaColor, 
		taxaNames = tree$tip.label,
		transparency = transparency
		)	
	##############################################
	# plot a tree but with blank tip labels
		# set x.lim so plot x limits is * (1 + extraMargin)
		# where 1 is the tree height (effectively)
	# calculate new x.lim by
		# *not* plotting a tree
	outPlot <- plot.phylo(
		tree,
		plot=FALSE,
		show.tip.label=FALSE,
		no.margin = removeSurroundingMargin,
		...)
	old_xlim <- outPlot$x.lim[2]
	new_xlim <- old_xlim * (1 + extraMargin)		
	#####
	par(new = TRUE)
	#####
	plot.phylo(
		tree,
		x.lim = new_xlim,
		show.tip.label = FALSE,
		no.margin = removeSurroundingMargin,
		...
		)
	##########################################
	# now get the last plotting environment
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	# get the device's aspect ratio
	devAspRatio <- grDevices::dev.size()[1] / grDevices::dev.size()[2]
	# get the plot's own aspect ratio
	plotAspRatio <- diff(lastPP$x.lim) / diff(lastPP$y.lim)
	# true aspect ratio is their product apparently
	plotAspRatio <- plotAspRatio / devAspRatio 
	##################################################
	#
	# pause 3 seconds so we don't spam the API
	#Sys.sleep(3)
	#	
	########################################
	for (i in 1:lastPP$Ntip){
		##########################
		# GET IMAGE
		picPNG <- getPhyloPicPNG(
			picID_PBDB = phylopicIDsPBDB[i], 
			cacheDir = cacheDir,
			cacheImage = cacheImage
			)
		######################################
		# PREP IMAGE
		picPNG <- prepPhyloPic(picPNG, 
			noiseThreshold = noiseThreshold,
			rescalePNG = rescalePNG, 
			trimPNG = trimPNG,
			makeMonochrome = makeMonochrome,
			plotComparison = FALSE
			)
		##################################
		# PLOT IMAGE
		plotSinglePhyloPic(
			picPNG = picPNG,
			whichTip = i,
			lastPP = lastPP,
			plotAspRatio = plotAspRatio,
			sizeScale = sizeScale,
			taxonColor = taxaColor[i]
			)	
		}
	modPhyloPlotInfo <- lastPP
	# add stuff here about what we did to the phylo plot?
		# like what?
	return(invisible(modPhyloPlotInfo))
	}


getPhyloPicPNG<-function(
		picID_PBDB, 
		cacheDir = NULL,
		cacheImage = FALSE
		){
	####################################################
	# first try to find and load a cached version
	# if that doesn't work
		# try to load from phylopic using PBDB UID
	# if that doesn't work
		# try to load the image from PBDB		
	notCached <- FALSE
	picPNG <- NULL
	###########################################
	# First try to get a cached version
	if(!is.null(cacheDir)) {
		cachePath <- file.path(cacheDir,
					paste0(picID_PBDB, ".png")
					)
		if(file.exists(cachePath)){
			picPNG <- png::readPNG(cachePath)
		}else{
			notCached <- TRUE
			}
		}
	##################################################
	# if that doesn't work
		# try to load from phylopic using PBDB UID
	if(is.null(picPNG)){
		picUID <- paste0("https://paleobiodb.org/data1.2/taxa/thumb.txt?id="
							,picID_PBDB)
		picUID <- read.csv(picUID, stringsAsFactors = FALSE)
		picUID <- picUID$uid
		picPNG <- getPhyloPicFromPhyloPic(picUID)
		}
	################################################
	# if that doesn't work
		# try to load the image from PBDB		
	if(is.null(picPNG)){
		picPNG <- getPhyloPicPNG_PBDB(picID_PBDB = picID_PBDB)
		}
	#########################
	if(cacheImage & notCached){
		if(!dir.exists(cacheDir)){
			dir.create(
				cacheDir,
				showWarnings = FALSE
				)
			}
		#
		png::writePNG(picPNG,
			target=file.path(
				cacheDir,
				paste0(picID_PBDB, ".png")
				)
			)
		}
	#########
	return(picPNG)	
	}		

	
plotSinglePhyloPic <- function(
		picPNG,
		whichTip,
		lastPP,
		sizeScale = 0.9,
		plotAspRatio,
		taxonColor = NULL
		){
		#
	#########################################
	#get aspect ratio
		# ratio of # of pixel dimensions
	picAspRatio <- dim(picPNG)[1]/dim(picPNG)[2]
	#############################################
	# adjustment of sizeScale
		# need to modify sizeScale relative to aspect ratio
	if(picAspRatio < 1){
		# its flatish so correct it by aspect ratio
		picSize <- sizeScale * (picAspRatio^0.7)
	}else{
		# then its skinny, not flat
		# don't do anything
		picSize <- sizeScale 
		}
	###########################################
	# GET THE COORDINATES
	#
	# offset is sizeScale/2 by default
	offset <- sizeScale*0.9* plotAspRatio
	#
	#points(lastPP$xx,lastPP$yy)	
	x<-lastPP$xx[whichTip]+offset
	y<-lastPP$yy[whichTip]
	#
	##################################################
	#adjust the position of the sides for the image
	#
	xAdj <- (picSize/2) * (plotAspRatio/picAspRatio) 
	yAdj <- picSize /2
	#
	##################################################
	# plot the picPNG using graphics::rasterImage
	picPNG_raster <- grDevices::as.raster(picPNG)
	# want color?
		# replaces all black tiles with the color in taxonColor
	if(!is.null(taxonColor)){
		print(taxonColor)
		picPNG_raster[which(picPNG_raster=="#000000FF")] <- taxonColor
		}
	# now plot the phylopic
	graphics::rasterImage(picPNG,	
		xleft = x - xAdj ,
		ybottom = y - yAdj ,
		xright = x + xAdj ,
		ytop = y + yAdj,
		interpolate = TRUE
		)
	# cool	
	}	

prepPhyloPic<-function(
		picPNG, 
		noiseThreshold = 0.1,
		rescalePNG = TRUE, 
		trimPNG = TRUE,
		makeMonochrome = FALSE,
		plotComparison = FALSE){
	############################################
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