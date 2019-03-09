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
#' placed around the plotted tree, or not. By default, \code{plotPhyloPicTree} will
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

# @param makeMonochrome If \code{TRUE}, PhyloPic silhouettes are
# forced to be purely monochrome black-and-white, with no gray
# scale. Most of the silhouettes are binary black-and-white already
# but some aren't, but those gray-scale values (sometimes?) seem
# to exist to indicate very fine features. However, maybe an image
# is far too much gray-scale, in which case users can try this
# option to force all silhouettes to be monochrome.
# The default is \code{FALSE}.

#' @param colorGradient Controls the depth gradient of color for the PhyloPics.
#' For typical plotting in black color, this means adjusting
#' the grayscale (and possibly removing any gray scale). 
#' Most of the silhouettes are binary black-and-white already but some
#' aren't, but those gray-scale values (sometimes?) seem
#' to exist to indicate very fine features. However, maybe an image
#' is far too much gray-scale, in which case users can apply this argument.
#' If \code{colorGradient = NULL} (the default), then nothing is adjusted.
#' If \code{colorGradient = "trueMonochrome"}, the entire image's gradients are
#' simplified to a duality: either fully colored or fully transparent.
#' If \code{colorGradient = "increaseDisparity"}, then a slightly less
#' extreme option is applied, with values transformed to greatly remove
#' in-between gray-scale value, shifting them toward color or
#' not-color without making the sihoullete purely monochrome.

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
#' individually. All PhyloPics expressed in colors other than the default
#' black are transformed as under the argument \code{colorGradient = "trueMonochrome"},
#' so that the PhyloPic is expressed with no intermediate gray-scale values.

#' @param orientation Controls the direction the phylogeny is plotted
#' in - can be either "rightwards" or "upwards".

#' @param transparency A numeric value between 0 and 1, either length 1, or the same
#' length as the number of tips on \code{tree}. This indicates the transparency of
#' either all the plotted PhyloPics, or allows user control over each PhyloPic
#' individually. The default is 1, which represents maximum opaqueness,
#' applied to all PhyloPics.
	
#' @param ... Additional arguments, passed to
#' \code{plot.phylo} for plotting of the tree. These
#' additional arguments may be passed to \code{plot},
#' and from there to \code{plot}. 
#' Some arguments are reserved and cannot be passed,
#' particularly: \code{direction}, \code{show.tip.label},
#' \code{no.margin}, \code{plot}, \code{xlim}, and\code{ylim}.

#' @return
#' This function silently returns the positions for elements in the tree

#' @seealso
#' See \code{\link{getTaxaDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhyloPicTree}}.

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
#' plotPhyloPicTree(tree = tree, 
#' 	taxaDataPBDB = taxaData)
#' 
#' }
#' 

# note sizeScale is vertical, proportional to the space between tips
	# max horizontal sizeScale stops flat / long PhyloPics from becoming overly huge

# set x.lim so plot x limits is * (1 + extraMargin)
# where 1 is the tree height (effectively)




#' @name plotPhyloPicTree
#' @rdname plotPhyloPicTree
#' @export
plotPhyloPicTree <- function(
		tree, 
		taxaDataPBDB = tree$taxaDataPBDB,
		# phylopicIDsPBDB = NULL, 
		#######################
		sizeScale = 0.9,
		removeSurroundingMargin = TRUE,
		extraMargin = 0.08,
		orientation = "rightwards",
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
		colorGradient = NULL,
		...
		){		
	#########################################
	# uses calls to the Paleobiology Database's API
		# or the phylopic API
		# to construct a phylogeny with PhyloPics
	# used as pictorial replacements for the tip labels
	# images taken directory from Phylopic, or PBDB
	###############################################
	# check that reserved arguments are not in ...
	reservedPlotArgs <- c("direction", "show.tip.label",
		"no.margin", "plot", "xlim", "ylim")
	dotArgNames <- names(list(...))
	if(any(!is.na(match(reservedPlotArgs, dotArgNames)))){
		matchingReserved <- match(reservedPlotArgs, dotArgNames)[
			!is.na(match(reservedPlotArgs, dotArgNames))
			]
		matchingReserved <- dotArgNames[matchingReserved]
		stop(paste0(
			"arguments '",
				paste0(matchingReserved, collapse= "' & '"),
				"' are reserved and cannot be directly passed to plot.phylo; see documentation"
			))
		}
	#############################################
	# test that orientation is one of the two values it can be
	if(length(orientation)!=1){
		stop("argument orientation must be a single value")}
	if(orientation != "upwards" & orientation != "rightwards"){
		stop('only orientation values of "upwards" and "rightwards" are currently accepted')
		}
	#################################################
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
	# get the device's aspect ratio
	devAspRatio <- grDevices::dev.size()[1] / grDevices::dev.size()[2]
	#
	# plot a tree but with blank tip labels
		# set x.lim so plot x limits is * (1 + extraMargin)
		# where 1 is the tree height (effectively)
	# calculate new x.lim by
		# *not* plotting a tree
	outPlot <- plot.phylo(
		tree,
		direction = orientation,
		plot=FALSE,
		show.tip.label=FALSE,
		no.margin = removeSurroundingMargin,
		...)
	# modify margins based on orientation
	if(orientation == "upwards"){
		# adjust extraMargin by aspect ratio
		extraMargin <- extraMargin/(devAspRatio^0.7)
		#
		new_xlim <- outPlot$x.lim[2] 
		new_ylim <- outPlot$y.lim[2] * (1 + extraMargin)
		}
	if(orientation == "rightwards"){
		# adjust extraMargin by aspect ratio
		extraMargin <- extraMargin*(devAspRatio^0.7)
		#
		new_xlim <- outPlot$x.lim[2] * (1 + extraMargin)
		new_ylim <- outPlot$y.lim[2] 
		}
	#####
	par(new = TRUE)
	#####
	plot.phylo(
		tree,
		x.lim = new_xlim,
		y.lim = new_ylim,
		direction = orientation,
		show.tip.label = FALSE,
		no.margin = removeSurroundingMargin,
		...
		)
	#
	##########################################
	# now get the last plotting environment
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	#
	# get the plot's own aspect ratio
	plotAspRatio <- diff(lastPP$x.lim) / diff(lastPP$y.lim)
	#
	# true aspect ratio is their product apparently
	plotAspRatio <- plotAspRatio / devAspRatio 
	#
	# calculate offset as a function of extraMargin
	offset <- (1+(extraMargin*0.5))
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
		# if this pic is colored, make it truly monochrome
		if(is.na(taxaColor[i])){
			colorGradientTaxon <- colorGradient
		}else{
			colorGradientTaxon <- "trueMonochrome"
			}
		#
		picPNG <- prepPhyloPic(picPNG, 
			noiseThreshold = noiseThreshold,
			rescalePNG = rescalePNG, 
			trimPNG = trimPNG,
			colorGradient = colorGradientTaxon,
			plotComparison = FALSE
			)
		##################################
		# PLOT IMAGE
		plotSinglePhyloPic(
			picPNG = picPNG,
			whichTip = i,
			lastPP = lastPP,
			offset = offset,
			orientation = orientation,
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
