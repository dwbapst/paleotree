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

#' @param sizeScale The default is \code{sizeScale = 0.9}.

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


# @param phylopicIDsPBDB ID numbers for images from Phylopic,
# as given by the Paleobiology Database's API under the output
# \code{image_no} (given when \code{show = img}).

#' @param ... Additional arguments, passed to
#' \code{plot.phylo} for plotting of the tree. These
#' additional arguments may be passed to \code{plot},
#' and from there to \code{plot}.


#' @param cacheDir If not \code{NULL}, first look here for a cached
#' version of the images. This makes loading faster.
#' The default is \code{NULL}.

# @param taxaColor If not \code{NULL},

# if NULL, all are black
# type character
# if its length 1
# if its a value that matches a tip label, color that taxon "red"
# if its a value that does not match a tip label, coerce that colors
	# if not colors, FAIL
# if its not length 1, it must be same length as number of tips
	# if not same length as number of tips, FAIL
# in which case each value color is expected
	# if not a color, FAIL
	

#### would be better to treat this as a vector of same length as
#### number of tip taxa, for which to indicate colors of

# @param ... Other arguments to pass to plotting code



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
	# max horizontal sizeScale stops flat / long phylopics from becoming overly huge

# set x.lim so plot x limits is * (1 + extraMargin)
# where 1 is the tree height (effectively)

# noiseThreshold threshold for noise in the PNG from Phylopic to be
# treated as meaningless noise (i.e. a color that is effectively whitespace)
# and thus can be trimmed as margin to be trimmed by the function


#' @name plotPhylopicTreePBDB
#' @rdname plotPhylopicTreePBDB
#' @export
plotPhylopicTreePBDB <- function(
		tree, 
		taxaDataPBDB = tree$taxaDataPBDB,
		# phylopicIDsPBDB = NULL, 
		sizeScale = 0.9,
		noiseThreshold = 0.1,
		extraMargin = 0.2,
		rescalePNG = TRUE,
		trimPNG = TRUE,
		makeMonochrome = FALSE,
		cacheDir = NULL,
		taxaColor = NULL,
		...
		){		
	#########################################
	# uses calls to the Paleobiology Database's API
		# to construct a phylogeny with phylopics
			# (from Phylopic, duh)
		# as pictorial replacements for the tip labels
	###############################################
	# check or obtain the phylopic IDs from PBDB
	phylopicIDsPBDB <- getPhyloPicIDNumFromPBDB(
		taxaData = taxaDataPBDB,
		#phylopicIDsPBDB = phylopicIDsPBDB,
		tree = tree)
	###############################################
	# determine colors for every taxon using taxaColor
	
		
	# if taxaColor is NULL, all are 'black'
		# just make it black
		# make all taxonColor NULL
		taxaColor <- rep(NULL,Ntip(tree))
	
	# else, taxaColor must be type character
		# if not, FAIL
	# if its length 1
		# if its a value that matches a tip label, color that taxon "red"
	
		#this code will make a single 'focal' taxon a bright red	
			# red is "#FF0000FF"
		
		
		# if its a value that does not match a tip label, coerce to a color
			# if not colors, FAIL
			# if colors, all taxa will be in that color
			# coerce to a hex value
			
	# if its not length 1, it must be same length as number of tips
		# if not same length as number of tips, FAIL
		# if right length, value is expected to be a color
			# if not a color, FAIL
		# if colors, each taxa will be in that color, in same order as tip.labels
			# coerce to a hex value		

	
	
	
	##############################################
	# plot a tree but with blank tip labels
		# set x.lim so plot x limits is * (1 + extraMargin)
		# where 1 is the tree height (effectively)
	# calculate new x.lim by
		# *not* plotting a tree
	outPlot <- 	plot.phylo(
		tree,
		x.lim = new_xlim,
		show.tip.label = FALSE,
		...
		)
	old_xlim <- outPlot$x.lim[2]
	new_xlim <- old_xlim * (1 + extraMargin)
	#####
	par(new = TRUE)
	#####
	plot.phylo(
		tree,
		x.lim = new_xlim,
		show.tip.label = FALSE,
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
	#
	
	
	
	# pause 3 seconds so we don't spam the API
	#Sys.sleep(3)
	
	
	
	for (i in 1:lastPP$Ntip){
		##########################
		# GET IMAGE
		picPNG <- getPhyloPicPNG(
			picID_PBDB = phylopicIDsPBDB[i], 
			cacheDir = cacheDir
			)
		######################################
		# PREP IMAGE
		picPNG <- prepPhyloPic(picPNG, 
			noiseThreshold = noiseThreshold,
			rescalePNG = rescalePNG, 
			trimPNG = trimPNG,
			makeMonochrome = makeMonochrome,
			plotComparison = plotComparison)
			}
		##################################
		# PLOT IMAGE
		plotSinglePhyloPic(
			picPNG = picPNG,
			whichTip = i,
			lastPP = lastPP,
			sizeScale = sizeScale,
			taxonColor = taxaColor[i]
			)	
				
		
		

		}
	modPhyloPlotInfo <- lastPP
	# add stuff here about what we did to the phylo plot
	


	
	return(invisible(modPhyloPlotInfo))
	}



getPhyloPicPNG<-function(
		picID_pbdb, 
		cacheDir = NULL
		){
	# first try to find and load a cached version
	# if that doesn't work
		# try to load from phylopic using PBDB UID
	# if that doesn't work
		# try to load the image from PBDB		
	cacheLater <- FALSE
	picPNG <- NULL
	###########################################
	# First try to get a cached version
	if(!is.null(cacheDir)) {
		cachePath <- file.path(cacheDir,
					paste0(picID, ".png")
					)
		if(file.exists(cachePath)){
			picPNG <- png:readPNG(cachePath)
		}else{
			cacheImage <- TRUE
			}
		}
	##################################################
	# if that doesn't work
		# try to load from phylopic using PBDB UID
	if(is.null(picPNG)){
		picUIDdataTable <- getPhyloPicUIDdataFromPBDB(picID = picID_pbdb)	
		picUID <- picUIDdataTable$uid
		picPNG <- getPhyloPicFromPhyloPic(picUID)
		}
	################################################
	# if that doesn't work
		# try to load the image from PBDB		
	if(is.null(picPNG)){
		picPNG <- getPhyloPicPNG_PBDB(picID_PBDB = picID_pbdb)
		}
	#########################
	if(cacheImage){
		
		}
		
	#########
	return(picPNG)	
	}		



	
	
plotSinglePhyloPic <- function(
		picPNG,
		whichTip,
		lastPP,
		sizeScale = 0.9,
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



	



	