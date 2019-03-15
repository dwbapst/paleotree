#' Date a Taxon-Topology from the Paleobiology Database Using Appearance Data from the API
#' 
#' The function \code{dateTaxonTreePBDB} takes a input consisting of
#' a topology, with tip and internal node labels corresponding to
#' taxa in the Paleobiology Database, and a table of data (containing those same tip and
#' node taxa) obtained from the taxa-list functionality of the Paleobiology Database's API,
#' with appearance times. This function will then output a tree with nodes reflecting the
#' ages of the respective higher taxa, based on their earliest times of appearance
#' from the Paleobiology Database.

#' @details
#' The dating by 







#' @param taxaTree A tree with tip taxon names matching the taxon names
#' in \code{taxaDataPBDB}. Probably a taxon tree estimated
#' with \code{\link{makePBDBtaxonTree}}.

#' @param taxaDataPBDB A data table of taxonomic information obtained
#' using the Paleobiology Database's API for a set of taxa that
#' includes the tip taxa on \code{taxaTree}, generated with
#' parameter \code{show = app} so that appearance times are included.

#' @param minBranchLen Following dating using the appearance times taken directly
#' from the PBDB for each tip taxon and node, the tree may then be assessed with
#' the minimum branch length algorithm, as applied by \code{\link{minBranchLength}}.
#' If \code{minBranchLen = 0}, the default, this step is skipped. It may be necessary
#' to set \code{minBranchLen} higher than zero to differentiate nodes in cases
#' with poor stratigraphic congruency, so that derived taxa are the first taxa
#' observed in a group.

#' @param plotTree If \code{TRUE}, the resulting dated tree is plotted.
#' This is \code{FALSE} by default.

#' @param tipTaxonDateUsed Controls what date for a taxon from the PBDB
#' is used for 'when' the tip should be placed in the dated phylogeny
#' produced by this function. The default, \code{tipTaxonDateUsed = "shallowestAge"}
#' will use the minimum age of the last appearance time of that taxon, which if it
#' is extant will be 0, and if it is extinct, will be the maximum constraint on the
#' age of its last appearance (i.e. the last time we saw it before it went extinct).
#' Other options are "deepestAge", which is the oldest possible first appearance time
#' from the PBDB, i.e. the maximum age constraint for the first appearance. As closely
#' related taxa often first occur in the same short interval of geologic time, due to
#' diversification bursts and/or the heterogeneity of fossil preservation, this may
#' result in confusing polytomies of many terminal taxa with no terminal branch lengths.

#' @return
#' Returns a dated phylogeny of class \code{phylo}, with an additional element
#' \code{$taxaDataPBDB} added containing the input \code{taxaDataPBDB}, as this might be
#' called by other functions.

#' @seealso
#' See \code{\link{getTaxaDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhyloPicTree}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.

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
#' tree <- makePBDBtaxonTree(data, rank = "genus") 
#' 
#' #get the ranges 
#' timeTree <- dateTaxonTreePBDB(tree)
#' 
#' plotPhyloPicTree(tree = timeTree)
#' axisPhylo()
#' 
#' }
#' 
#' \dontrun{
#' 
#' # example using strap
#' 
#' rangesFourDate<-timeTree$tipTaxonFourDateRanges
#' colnames(rangesMinMax )<-c("FAD","LAD")
#' 
#' library(strap)
#' #now plot it
#' strap::geoscalePhylo(
#' 	tree=timeTree,
#' 	direction = "upwards",
#' 	ages=rangesMinMax,
#' 	cex.tip=0.7,
#' 	cex.ts=0.55,
#' 	cex.age=0.5,
#' 	width=3,
#' 	tick.scale = 50,
#' 	quat.rm=TRUE,
#' 	boxes = "Period",
#' 	arotate = 90,
#' 	units=c("Eon","Period","Era"),
#' 	x.lim=c(650,-20))
#' 
#' }
#' 


#' @name dateTaxonTreePBDB
#' @rdname dateTaxonTreePBDB
#' @export
dateTaxonTreePBDB <- function(
		taxaTree,
		taxaDataPBDB = taxaTree$taxaDataPBDB,
		minBranchLen = 0,
		tipTaxonDateUsed = "shallowestAge",
		plotTree = FALSE){
	###################################
	if(!any(colnames(taxaDataPBDB)!="lastapp_min_ma")){
		stop(paste0(
			"input taxaTree must have a taxaDataPBDB element of PBDB data",
			"generated with show=app"))
		}
	############################
	# replace min & max last appearance ages
		# with 0 if "is_extant" is "extant"
	isExtant <- taxaDataPBDB$is_extant == "extant"
	lastAppTimes <- c("lastapp_min_ma","lastapp_max_ma")
	firstAppTimes <- c("firstapp_min_ma","firstapp_max_ma")
	taxaDataPBDB[isExtant, lastAppTimes] <- 0
	# if both first appearance times for a taxon are NA (no fossil occurrences)
		# and its extant, put in 0 for its first appearance time
	hasNAFirstApps <- apply(taxaDataPBDB[,firstAppTimes],1,
		function(x) is.na(x[1]) & is.na(x[2])
		)
	taxaDataPBDB[isExtant & hasNAFirstApps, firstAppTimes] <- 0
	###########################################
	# get four date taxon ranges for all tip taxa
	# match tip-taxa to taxa-data
		#sort based on tip labels
	tipMatches <- match(taxaTree$tip.label, taxaDataPBDB$taxon_name)
	#
	tipTaxonFourDateRanges <- taxaDataPBDB[tipMatches,
		c("firstapp_max_ma","firstapp_min_ma",
			"lastapp_max_ma","lastapp_min_ma")]
	# select the right tipAges based on tipTaxonDateUsed
	if(tipTaxonDateUsed == "shallowestAge"){
		tipAges <- tipTaxonFourDateRanges$lastapp_min_ma
		}
	if(tipTaxonDateUsed == "deepestAge"){
		tipAges <- tipTaxonFourDateRanges$firstapp_max_ma
		}
	#
	#############################################
	# get node ages
	nodeNames<-paste0(taxaTree$node.label,collapse=",")
	apiAddressNodes <- paste0(
		"http://paleobiodb.org/data1.2/taxa/list.txt?name=",
		nodeNames,"&show=app,parent"
		)
	# browseURL(apiAddressNodes)
	nodeData <- read.csv(apiAddressNodes ,
		stringsAsFactors = FALSE)
	# get node max ages
	nodeMaxAges <- nodeData$firstapp_max_ma[
		match(taxaTree$node.label, nodeData$taxon_name)]
	# construct tip+node age vector
	allAges <- as.numeric(c(tipAges, nodeMaxAges))
	##########################################
	# get dated tree
	datedTree <- nodeDates2branchLengths(
		nodeDates = allAges,
		tree = taxaTree,
		allTipsModern = FALSE)
	#
	############################################
	# check that the tree and its root age makes sense
	#checkRes <- checkRootTime(datedTree)			
	#################################
	if(minBranchLen>0){
		# take care of zero length branches
		datedTree<-minBranchLength(
			tree = datedTree,
			mbl = minBranchLen)
		datedTree<-ladderize(datedTree)
		#
		plotName <- paste0("Dated Phylogeny (",
			minBranchLen,
			" Minimum Branch Length)"
			)
	}else{
		plotName <- "Dated Phylogeny"
		}
	#
	############################################
	# check that the tree and its root age makes sense
	checkRootTimeRes <- checkRootTime(datedTree,
		stopIfFail = TRUE)
	#################
	if(plotTree){
		plot(datedTree, main=plotName,
			show.node.label=FALSE, cex=0.5)
		axisPhylo()
		}
	###################		
	# output four date taxon ranges for all tip taxa
		# (for strap or otherwise)	
	rownames(tipTaxonFourDateRanges) <- taxaTree$tip.label
	datedTree$tipTaxonFourDateRanges <- tipTaxonFourDateRanges
	# output the original taxaDataPBDB
	datedTree$taxaDataPBDB <- taxaDataPBDB
	# return
	return(datedTree)
	}