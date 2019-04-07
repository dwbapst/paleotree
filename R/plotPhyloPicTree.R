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
#' construct a taxonomy-tree. If the functions listed in \code{\link{getDataPBDB}}
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

#' @param maxAgeDepth The maximum tree depth displayed for a tree given with
#' branch lengths (age depth for a dated tree). The portion of the phylogeny
#' older than this date will not be shown. \code{NULL} by default. If provided,
#' the input tree must have branch lengths in \code{tree$edge.length}.

#' @param noiseThreshold A threshold for noise in the PNG from PhyloPic
#' to be treated as meaningless noise (i.e. a color that is effectively
#' whitespace) and thus can be trimmed as empty margin which can be
#' trimmed before the silhouette is plotted. The units for this argument
#' are on a scale from 0 to 1, with 0 being true white space, and values
#' between 0 and 0.5 representing colors closer to whitespace than true
#' black. The default is \code{noiseThreshold = 0.1}.

# noiseThreshold threshold for noise in the PNG from PhyloPic to be
# treated as meaningless noise (i.e. a color that is effectively whitespace)
# and thus can be trimmed as margin to be trimmed by the function

#' @param removeSurroundingMargin This argument controls the \code{no.margin} argument
#' in the function \code{plot.phylo}, which controls whether a (very large) margin is 
#' placed around the plotted tree, or not. By default, \code{plotPhyloPicTree} will
#' suppress that margin, so that the plotted tree goes (very nearly) to the edges
#' of the plotting area.

# @param extraMargin How much extra margin should be added to the side of the graph
# to which PhyloPics are being added? This value is 0.08 by default, which works well
# if the margin surrounding the entire plot is suppressed via argument
# \code{removeSurroundingMargin = TRUE}. If the surrounding margin is not suppressed,
# plots look mostly okay if users change to \code{extraMargin = 2}. Obviously this value
# should be tweaked for every tree, size of plot and aspect ratio for maximum clarity. 

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
#' If \code{colorGradient = "original"} (the default), then nothing is adjusted.
#' If \code{colorGradient = "trueMonochrome"}, the entire image's gradients are
#' simplified to a duality: either fully colored or fully transparent.
#' If \code{colorGradient = "increaseDisparity"}, then a slightly less
#' extreme option is applied, with values transformed to greatly remove
#' in-between gray-scale value, shifting them toward color or
#' not-color without making the silhouette purely monochrome.

# @param phylopicIDsPBDB ID numbers for images from PhyloPic,
# as given by the Paleobiology Database's API under the output
# \code{image_no} (given when \code{show = img}).

#' @param addTaxonStratDurations If \code{TRUE}, solid color boxes are
#' plotted on the tree to indicated known taxon ranges, from the oldest 
#' possible for the oldest known observation of that taxon, to the youngest 
#' possible age for the youngest known observation of that taxon. This 
#' information needs to be supplied as input, see argument \code{taxaStratRanges}.
#' If \code{FALSE} (the default), nothing happens.

#' @param stratDurationBoxWidth The width of the stratigraphic duration boxes
#' plotted for taxa on the tree. By default, this is 0.7 units. If
#' \code{addTaxonStratDurations = FALSE} (the default), this argument affects nothing.

#' @param taxaStratRanges A matrix of four-date range information, as is often
#' used when converting Paleobiology Database taxon data to a dated tree. By
#' default, this is expected to be located at \code{tree$tipTaxonFourDateRanges},
#' which is where such data is placed by default by the function
#' \code{\link{dateTaxonTreePBDB}}. If \code{addTaxonStratDurations = FALSE}
#' (the default), this data is not checked for.

#' @param colorAxisPhylo A color in which the axis for the phylogenetic's depth
#' (generally a time-scale) will be plotted in, for both the axis, its
#' tickmarks, and the labels for the tickmarks.

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

#' @param depthAxisPhylo If \code{TRUE}, the \code{ape} function
#' \code{axisPhylo} is run to add an axis of the tree depth to
#' the tree, which must have branch lengths. \code{FALSE} by default.
#' If \code{removeSurroundingMargin = TRUE}, which removes extraneous margins,
#' the margins of the plot will be adjusted to make room for the plotted axis.

#' @param resetGrPar If \code{TRUE} (the default), the graphic parameters are
#' reset, so that choices of margins and coordinate system manipulation done
#' as part of this function do not impact the next plot made in this graphic
#' device. If you need to add additional elements to the plot after running
#' this function within R, you should set this argument to \code{FALSE}.

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
#' This function silently returns the positions for elements in the
#' tree (.e. the environmental information obtained about the
#' previous plotting environment of the tree as plotted), along
#' with a saved set of the graphic parameters as they were
#' at the end of the function's run.

#' @seealso
#' See \code{\link{getDataPBDB}}, \code{\link{makePBDBtaxonTree}},
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
#' library(paleotree)
#' 
#' taxaAnimals<-c("Archaeopteryx", "Eldredgeops",
#' 	"Corvus", "Acropora", "Velociraptor", "Gorilla", 
#' 	"Olenellus", "Lingula", "Dunkleosteus",
#' 	"Tyrannosaurus", "Triceratops", "Giraffa",
#' 	"Megatheriidae", "Aedes", "Histiodella",
#' 	"Rhynchotrema", "Pecten", "Homo", "Dimetrodon",
#' 	"Nemagraptus", "Panthera", "Anomalocaris")
#' 
#' data <-getSpecificTaxaPBDB(taxaAnimals)
#' tree <- makePBDBtaxonTree(data, rankTaxon = "genus") 
#' 
#' plotPhyloPicTree(tree = tree)
#' 
#' # let's plot upwards but at a funny size
#' dev.new(height = 5, width = 10)
#' plotPhyloPicTree(tree = tree,
#' 	 orientation = "upwards") 
#' 
#' # dated tree plotting
#' 
#' #date the tree 
#' timeTree <- dateTaxonTreePBDB(tree, minBranchLen = 10)
#' 
#' plotPhyloPicTree(tree = timeTree)
#' 
#' # plotting the dated tree with an axis
#' plotPhyloPicTree(tree = timeTree,
#' 	depthAxisPhylo= TRUE)
#' 
#' # now upwards!
#' plotPhyloPicTree(tree = timeTree,
#' 	orientation = "upwards",
#' 	depthAxisPhylo= TRUE)
#' 
#' ###################################
#' # plotting a time tree with stratigraphic ranges
#' 
#' plotPhyloPicTree(tree = timeTree,
#'    addTaxonStratDurations = TRUE)
#' 
#' plotPhyloPicTree(tree = timeTree,
#'    addTaxonStratDurations = TRUE,
#'    orientation = "upwards",
#'    depthAxisPhylo= TRUE)
#' 
#' ########
#' # adjusting a tree to ignore a very old root
#' 
#' # let's pretend that metazoans are extremely old
#' treeOldRoot <- timeTree
#' rootEdges <- timeTree$edge[,1] == (Ntip(timeTree)+1)
#' rootEdgeLen <- timeTree$edge.length[rootEdges]
#' treeOldRoot$edge.length[rootEdges] <- rootEdgeLen + 1500
#' treeOldRoot$root.time <- NULL
#' 
#' # plot it
#' plot(treeOldRoot)
#' axisPhylo()
#' # yep, that's really old
#' 
#' # let's plot it now with the PhyloPic
#' plotPhyloPicTree(tree = treeOldRoot,
#' 	depthAxisPhylo = TRUE)
#' 
#' # let's crop that old lineage
#' plotPhyloPicTree(tree = treeOldRoot,
#' 	maxAgeDepth = 500,
#' 	depthAxisPhylo = TRUE)
#' # cool!
#' 
#' ##################################
#' # playing with colors
#' plotPhyloPicTree(tree = tree,
#' 	 taxaColor = "green")
#' 
#' # inverting the colors
#' par(bg="black")
#' taxaColors <- rep("white",Ntip(tree))
#' # making a red giraffe
#' taxaColors[4] <- "red"
#' plotPhyloPicTree(
#' 	tree = tree, 
#' 	orientation = "upwards",
#' 	edge.color = "white",
#' 	taxaColor=taxaColors)
#' 
#' ######################################
#' # let's try some different phylopics
#' 	  # like a nice tree of commonly known tetrapods
#' 
#' tetrapodList<-c("Archaeopteryx", "Columba", "Ectopistes",
#' 	"Corvus", "Velociraptor", "Baryonyx", "Bufo",
#' 	"Rhamphorhynchus", "Quetzalcoatlus", "Natator",
#' 	"Tyrannosaurus", "Triceratops", "Gavialis",
#' 	"Brachiosaurus", "Pteranodon", "Crocodylus",
#' 	"Alligator", "Giraffa", "Felis", "Ambystoma",
#'  	"Homo", "Dimetrodon", "Coleonyx", "Equus",
#' 	"Sphenodon", "Amblyrhynchus")
#' 
#' data <-getSpecificTaxaPBDB(tetrapodList)
#' 
#' tree <- makePBDBtaxonTree(data, rankTaxon = "genus")
#' 
#' plotPhyloPicTree(tree = tree)
#' 
#' 
#' }
#' ####################################
#' \dontrun{
#' # let's check out speed increase from caching!
#'     # can try this on your own machine
#' 
#' #first time
#' system.time(plotPhyloPicTree(tree = tree))
#' # second time
#' system.time(plotPhyloPicTree(tree = tree))
#' 
#' ##################################
#' # make a pretty plot
#' 
#' taxaSeventyEight <- c(
#' 	"Archaeopteryx", "Pinus", "Procoptodon", "Olenellus", "Eldredgeops",
#' 	"Quetzalcoatlus", "Homo", "Tyrannosaurus", "Triceratops", "Giraffa",
#' 	"Bolivina", "Cancer", "Dicellograptus", "Dunkleosteus", "Solanum",
#' 	"Anomalocaris", "Climacograptus", "Halysites", "Cyrtograptus", 
#' 	"Procoptodon", "Megacerops", "Moropus", "Dimetrodon", "Lingula",
#' 	"Rhynchosaurus", "Equus", "Megaloceros", "Rhynchotrema", "Pecten",
#' 	"Echinaster", "Eocooksonia", "Neospirifer", # "Prototaxites", 
#' 	"Cincinnaticrinus", "Nemagraptus", "Monograptus", "Pongo", "Acropora",
#' 	"Histiodella", "Agathiceras", "Juramaia", "Opabinia", "Arandaspis",
#' 	"Corvus", "Plethodon", "Latimeria", "Phrynosoma", "Araucarioxylon",
#' 	"Velociraptor", "Hylonomus", "Elginerpeton", "Rhyniognatha",
#' 	"Tyto", "Dromaius", "Solenopsis", "Gorilla", "Ginkgo", "Terebratella", 
#' 	"Caretta", "Crocodylus", "Rosa", "Prunus", "Lycopodium", "Meganeura",
#' 	"Diplodocus", "Brachiosaurus", "Hepaticae", "Canadaspis", "Pikaia",
#' 	"Smilodon", "Mammuthus", "Exaeretodon", "Redondasaurus", "Dimetrodon",
#' 	"Megatheriidae", "Metasequoia", "Aedes", "Panthera", "Megalonyx")
#' 
#' data <-getSpecificTaxaPBDB(taxaSeventyEight)
#' tree <- makePBDBtaxonTree(data, rankTaxon = "genus") 
#' 
#' timeTree <- dateTaxonTreePBDB(tree,
#'   minBranchLen = 10)
#' 
#' date <- format(Sys.time(), "%m-%d-%y")
#' file <- paste0(
#'	"tree_taxa78_phylopic_stratTree_",
#' 	date, ".pdf")
#' 
#' png(file = file,
#' 	height = 5, width = 12, 
#' 	units = "in", res = 300)
#' par(bg="black")
#' par(mar=c(0,0,3,0))
#' taxaColors <- rep("white", Ntip(timeTree))
#' taxaColors[4] <- "red"
#' 
#' plotPhyloPicTree(
#' 	tree = timeTree, 
#' 	orientation = "upwards",
#' 	addTaxonStratDurations = TRUE,
#' 	edge.color = "white",
#' 	maxAgeDepth = 700,
#' 	taxaColor=taxaColors,
#' 	depthAxisPhylo = TRUE,
#' 	colorAxisPhylo = "white")
#' dev.off()
#' shell.exec(file)
#' 
#' }
#' 

#' @name plotPhyloPicTree
#' @rdname plotPhyloPicTree
#' @export
plotPhyloPicTree <- function(
		tree, 
		taxaDataPBDB = tree$taxaDataPBDB,
		# phylopicIDsPBDB = NULL, 
		#extraMargin = 0.08,
		#######################
		maxAgeDepth = NULL,
		depthAxisPhylo = FALSE,
		colorAxisPhylo  = "black",
		#######################
		addTaxonStratDurations = FALSE,
		taxaStratRanges = tree$tipTaxonFourDateRanges,
		stratDurationBoxWidth = 0.7,
		####################
		sizeScale = 0.9,
		removeSurroundingMargin = TRUE,
		orientation = "rightwards",
		resetGrPar = TRUE,
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
		colorGradient = "original",
		...
		){		
	#########################################
	# uses calls to the Paleobiology Database's API
		# or the phylopic API
		# to construct a phylogeny with PhyloPics
	# used as pictorial replacements for the tip labels
	# images taken directory from PhyloPic, or PBDB
	###############################################
	# save original graphic par
	oldPar<-par(no.readonly = TRUE)
	##########################
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
	############################################
	# check that the tree and its root age makes sense
	checkRes <- checkRootTime(tree)			
	#############################
	# get ranges if addTaxonStratDurations
		# then need to adjust tree so that plotting takes into account range ages
	if(addTaxonStratDurations){
		if(is.null(tree$edge.length)){
			stop("addTaxonStratDurations is inapplicable for undated trees - input tree has no branch lengths")
		}else{
			if(is.null(taxaStratRanges)){
				stop("addTaxonStratDurations = TRUE, but taxaStratRanges not supplied")
				}
			ranges <- getSortedMinMaxStratRanges(timeTree = tree,
				rangesFourDate = taxaStratRanges)
			}
	}else{
		ranges<-NULL	
		}
	#############################
	# get root age if not already present
	if(!is.null(tree$edge.length) & is.null(tree$root.time)){
		# not sure how this could happen but...
		tree$root.time <- max(node.depth.edgelength(tree))
		}
	#############################
	# get maximum age depth, which should be a function of maxAgeDepth and tree depth only
	if(is.null(maxAgeDepth)){
		if(is.null(tree$edge.length)){
			maxAgeLim <- NA
		}else{
			maxAgeLim <- tree$root.time
			}
	}else{
		if(length(maxAgeDepth) != 1 | !is.numeric(maxAgeDepth)){
			stop("maxAgeDepth must be a numeric vector of length 1 if provided")
			}
		# test if the tree even has edge lengths
			# if it does, see if the edge lengths mean anything 
		if(is.null(tree$edge.length)){
			stop("maxAgeDepth provided but input tree has no branch lengths, and thus is undated - please reconcile")
		}else{
			# in absolute time
			maxAgeLim <- maxAgeDepth
			}
		}
	###################################
	# now need to get minAgeLim
		# need to adjust tree so that plotting takes into account range ages	
	if(is.null(ranges)){
		if(is.null(tree$edge.length)){
			minAgeLim <- NA
		}else{
			furthestTipDist <- max(node.depth.edgelength(tree))
			minAgeLim <- tree$root.time - furthestTipDist
			}		
	}else{
		minAgeLim <- min(ranges)
		}
	############
	if(is.na(maxAgeLim) & is.na(minAgeLim)){
		ageLimInput <- NULL
	}else{
		ageLimInput <- c(maxAgeLim, minAgeLim)
		# but now need to make from root = 0
			# subtract from the age of the root, which will become zero		
		ageLimInput<- tree$root.time - c(maxAgeLim, minAgeLim)
		# rescale ranges relative to the maximum age depth
		rangesRescaled <- tree$root.time - ranges	
		}
	##################
	# assign limits depending on orientation
	if(orientation == "rightwards"){
		xlimInput <- ageLimInput
		ylimInput <- NULL
		}
	if(orientation == "upwards"){
		xlimInput <- NULL
		ylimInput <- ageLimInput		
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
	############################################
	if(depthAxisPhylo & is.null(tree$edge.length)){
		stop(paste0("depthAxisPhylo is TRUE but tree doesn't have branch lengths -",
			"\nWhy do you want to plot a phylo axis on a tree with no branch lengths??"
			))
		}
	if(removeSurroundingMargin & depthAxisPhylo){
		# reset removeSurroundingMargin
		removeSurroundingMargin <- FALSE
		if(orientation == "rightwards"){
			par(mar = c(2.5,0,0,0))
			}
		if(orientation == "upwards"){
			par(mar = c(0,0,0,2.5))		
			}		
		}
	##############################################
	# plot a tree but with blank tip labels
		# so can adjust x.lim to include an extra margin
	# calculate new x.lim by
		# *not* plotting a tree
	outPlot <- plot.phylo(
		tree,
		show.tip.label=FALSE,
		x.lim = xlimInput,
		y.lim = ylimInput,
		direction = orientation,
		no.margin = removeSurroundingMargin,
		plot=FALSE,
		...)
	#
	# get the device's aspect ratio
	devAspRatio <- grDevices::dev.size()[1] / grDevices::dev.size()[2]
	#
	# get the x and y lims of the previous plot
	new_xlim <- c(outPlot$x.lim[1], 
		outPlot$x.lim[2])	
	new_ylim <- c(outPlot$y.lim[1], 
		outPlot$y.lim[2])	
	plotSizeX <- new_xlim[2] - new_xlim[1] 
	plotSizeY <- new_ylim[2] - new_ylim[1]
	#
	# could also get from par("usr")
	#plotSizeX <- par("usr")[2] - par("usr")[1]
	#plotSizeY <- par("usr")[4] - par("usr")[3]
	#
	# get how many user coord units are in inches for each axis
	OneUserCoordUnitX_in_InchesX <- grDevices::dev.size()[1] / plotSizeX
	OneUserCoordUnitY_in_InchesY <- grDevices::dev.size()[2] / plotSizeY
	#
	# calculate how many horizontal user coords are equivalent to how many vertical user coord
	OneUserCoordUnitY_in_UserCoordUnitX <- OneUserCoordUnitY_in_InchesY / OneUserCoordUnitX_in_InchesX
	# calculate how many horizontal user coords are equivalent to how many vertical user coord
	OneUserCoordUnitX_in_UserCoordUnitY <- OneUserCoordUnitX_in_InchesX / OneUserCoordUnitY_in_InchesY
	#	
	# modify margins based on orientation
	if(orientation == "rightwards"){
		# extra margin is one vertical unit, in horizontal user coord units
			# adjust by aspect ratio and sizescale/2
		extraMargin <- OneUserCoordUnitY_in_UserCoordUnitX * (sizeScale/2) / devAspRatio
		new_xlim[2] <- new_xlim[2] + extraMargin
		}
	if(orientation == "upwards"){
		# extra margin is one horizontal unit, in vertical user coord units
			# adjust by aspect ratio and sizescale/2
		extraMargin <- OneUserCoordUnitX_in_UserCoordUnitY * (sizeScale/2) * devAspRatio
		new_ylim[2] <- new_ylim[2] + extraMargin
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
	# plot axisPhylo if depthAxisPhylo 
	if(depthAxisPhylo){
		if(orientation == "rightwards"){
			axisPhylo(side=1,
				col = colorAxisPhylo, col.axis = colorAxisPhylo)
			}
		if(orientation == "upwards"){
			axisPhylo(side=4,
				col = colorAxisPhylo, col.axis = colorAxisPhylo)
			}		
		}
	##########################################
	# now get the last plotting environment
	lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
	############################################
	if(addTaxonStratDurations){
		#
		newXY <- plotTaxonStratDurations(
			rangesMinMax = rangesRescaled,
			orientation = orientation,
			XX = lastPP$xx, YY = lastPP$yy,
			boxWidth = stratDurationBoxWidth, 
			boxCol = taxaColor)
	}else{
		newXY <- list(XX = lastPP$xx,
			YY = lastPP$yy)
		}
	#########################################################
	# get the plot's own aspect ratio
	#plotAspRatio <- diff(lastPP$x.lim) / diff(lastPP$y.lim)
	# actually... better to get from par("usr")
		# xmin, ymin, xmax, ymax
	plotDimensions <- list(
		xmin = par("usr")[1],
		xmax = par("usr")[2],
		ymin = par("usr")[3],
		ymax = par("usr")[4]
		)
	plotSizeX <- plotDimensions$xmax - plotDimensions$xmin
	plotSizeY <- plotDimensions$ymax - plotDimensions$ymin
	plotAspRatio <- plotSizeX / plotSizeY
	#
	# true aspect ratio is their product apparently
	plotAspRatio <- plotAspRatio / devAspRatio 
	#
	# calculate offsetPic as a function of orientation and plot limits
	if(orientation == "rightwards"){
		offsetPic <- plotDimensions$xmax - max(newXY$XX)
		}
	if(orientation == "upwards"){
		offsetPic <- plotDimensions$ymax - max(newXY$YY)
		}
	#
	offsetPic <- offsetPic * 0.5
	#
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
			XX = newXY$XX,
			YY = newXY$YY,
			offsetPic = offsetPic,
			orientation = orientation,
			plotAspRatio = plotAspRatio,
			sizeScale = sizeScale,
			taxonColor = taxaColor[i]
			)	
		}
	modPhyloPlotInfo <- lastPP
	# add stuff here about what we did to the phylo plot?
		# like what?
	# adjusted graphic parameters
	modPhyloPlotInfo$savedPar <- par(no.readonly = TRUE)
	# adjusted XY
	modPhyloPlotInfo$newXY <- newXY
	# reset graphic par
	if(resetGrPar){
		suppressWarnings(par(oldPar))
		}
	#
	return(invisible(modPhyloPlotInfo))
	}
