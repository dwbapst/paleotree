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
#' The dating by this function is very simplistic, representing a rather
#' straight interpretation of what the PBDB reports. The dated trees
#' produced should not be taken overly seriously.

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

#' @param dropZeroOccurrenceTaxa If \code{TRUE}, the default, then extinct taxa
#' or extinct clades found to have zero occurrences in the Paleobiology Database
#' are removed. If this option isn't used, the function will likely fail as nodes
#' or tips with \code{NA} ages listed cannot be processed by \code{parentChild2TaxonTree}.

#' @return
#' Returns a dated phylogeny of class \code{phylo}, with an additional element
#' \code{$taxaDataPBDB} added containing the input \code{taxaDataPBDB}, as this might be
#' called by other functions.

#' @seealso
#' See \code{\link{getDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhyloPicTree}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.
#' 
#' The equuid tree used in the examples is from:
#' MacFadden, B. J. 1992. Fossil horses: systematics, paleobiology, and
#' evolution of the family Equidae. \emph{Cambridge University Press}.

#' @examples
#' 
#' \donttest{
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
#' #get the ranges 
#' timeTree <- dateTaxonTreePBDB(tree)
#' 
#' plotPhyloPicTree(tree = timeTree,
#'      depthAxisPhylo = TRUE)
#' 
#' }
#' 
#' ####################################
#' 
#' \dontrun{
#' 
#' # can also plot dated tree with strap
#' 
#' library(strap)
#' #now plot it
#' strap::geoscalePhylo(
#'     tree=timeTree,
#'     direction = "upwards",
#'     ages=rangesMinMax,
#'     cex.tip=0.7,
#'     cex.ts=0.55,
#'     cex.age=0.5,
#'     width=3,
#'     tick.scale = 50,
#'     quat.rm=TRUE,
#'     boxes = "Period",
#'     arotate = 90,
#'     units=c("Eon","Period","Era"),
#'     x.lim=c(650,-20)
#'     )
#' }
#' 
#' ##############################################################
#' 
#' \donttest{
#' 
#' # we can also use this for pre-existing trees
#'     # for example, this tree of equuids (horses)
#'     # borrowed from UCMP materials on horse evolution
#'     # https://evolution.berkeley.edu/evolibrary/images/HorseTree.pdf
#'     # (apparently from MacFadden, 1992? Citation above)
#' 
#' # read the tree in as Newick string
#' horseTree <- ape::read.tree(file=NULL, 
#'     text = paste0(
#'          "(Eohippus,(Xenicohippus,(Haplohippus,(Epihippus,",
#'          "(Miohippus,(((Hypohippus,Megahippus),(Anchitherium,",
#'          "Kalobatippus)),(Archaeohippus,(Desmatippus,(Parahippus,",
#'          "(Merychippus,(((Hipparion_Merychippus,(Nannippus,",
#'          " Cormohipparion)),(Pseudhipparion,(Neohipparion,",
#'          " Hipparion))),(Equine_Merychippus,((Protohippus,Calippus),",
#'          "(Pliohippus,(Astrohippus,(Dinohippus,Equus))))))))))))))));"
#'          )
#'     )
#' 
#' # note there is a message that the tree lacks node names
#'     # this is unexpected / atypical for taxon trees
#' 
#' plot(horseTree)
#' 
#' # now let's get data on the tip from the PBDB
#'     # using getSpecificTaxaPBDB
#' horseData <- getSpecificTaxaPBDB(horseTree$tip.label)
#' 
#' # now we can date the tree with dateTaxonTreePBDB
#' 
#' datedHorseTree <- dateTaxonTreePBDB(
#'     taxaTree = horseTree,
#'     taxaDataPBDB = horseData,
#'     minBranchLen = 1
#'     )
#' 
#' # and let's try plotting it!	
#' plotPhyloPicTree(
#'     tree = datedHorseTree,
#'     depthAxisPhylo = TRUE
#'     )		
#' 	
#' # a fairly boring phylopic diagram
#'     # not many horse phylopics as of 07-16-19?
#' }
#' 
#' \dontrun{
#' 
#' # Let's look at this horse tree with strap
#'
#' library(strap)
#' 
#' geoscalePhylo(
#'     tree = datedHorseTree,
#'     ages = datedHorseTree$ranges.used,
#'     cex.tip = 0.7,
#'     cex.ts = 0.7,
#'     cex.age = 0.7,
#'     width = 4,
#'     tick.scale = 15,
#'     boxes = "Epoch",
#'     erotate = 90,
#'     quat.rm = TRUE,
#'     units = c("Period","Epoch"),
#'     x.lim = c(65,-10)
#'     )
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
		dropZeroOccurrenceTaxa = TRUE,
		plotTree = FALSE){
	###################################
	if(!any(colnames(taxaDataPBDB)!="lastapp_min_ma")){
		stop(paste0(
			"a data table of PBDB variables, generated with show=app",
			" must be provided either directly",
			"or as a taxaDataPBDB element of taxaTrees"))
		}
	############################
	lastAppTimes <- c("lastapp_min_ma","lastapp_max_ma")
	firstAppTimes <- c("firstapp_min_ma","firstapp_max_ma")
	#
	# save orig taxa data
	taxaDataPBDB_orig <- taxaDataPBDB
	#
	############################################################
	# 03-24-19
	# why don't I get all tip and node data simultaneously?
	#
	# skip if there aren't any node labels
	if(is.null(taxaTree$node.label)){
		#
		message("No node labels found - is this a taxon tree from the PBDB?")
		# just take original taxon data as the 'combined data'	
		combData <- taxaDataPBDB
		#
	}else{
		#
		# get node ages
		nodeNames <- taxaTree$node.label
		# remove any ".1" in the taxon names
			# hopefully there aren't any "." in the real taxon names... sigh...
		nodeNames <- sapply(
			strsplit(
				nodeNames, 
				split=".", 
				fixed=TRUE
				),
			function(x) x[[1]]
			)
		# remove nodeNames not in taxaDataPBDB
		nodeNamesNoMatch <- sapply(
			nodeNames, 
			function(x) 
				all(x != taxaDataPBDB$taxon_name)
				)
		if(any(nodeNamesNoMatch)){
			nodeNames <- nodeNames[nodeNamesNoMatch]
				# get API URL
			nodeNames <- paste0(nodeNames, collapse=",")
			apiAddressNodes <- paste0(
				"http://paleobiodb.org/data1.2/taxa/list.txt?name=",
				nodeNames,"&show=app,parent"
				)
			# browseURL(apiAddressNodes)
			nodeData <- read.csv(apiAddressNodes,
				stringsAsFactors = FALSE)
			# combine with taxon data
				# reducing scope to same columns as nodeData
					#print(colnames(nodeData))
					#print(colnames(taxaDataPBDB))
			sharedColNames <- intersect(colnames(nodeData), colnames(taxaDataPBDB))
			taxaDataReduced <-  taxaDataPBDB[,sharedColNames]
			nodeDataReduced <- nodeData[,sharedColNames]
			combData <- rbind(taxaDataReduced, nodeDataReduced)
			}
		}
	#		
	appData <- processTaxonAppData(dataPBDB = combData, 
		dropZeroOccurrenceTaxa = dropZeroOccurrenceTaxa)
	#
	#################
	# drop tips not in appData
	dropTips <- sapply(taxaTree$tip.label,
		function(x) all(x != appData$taxon_name)
		)
	# drop if any need to be dropped
	if(any(dropTips)){
		message(paste0("The following tips did not have resolvable",
			" appearance times and were dropped:"))
		message(paste0(taxaTree$tip.label[dropTips],
			collapse = ", "))
		taxaTree <- drop.tip(taxaTree,
			tip = which(dropTips)
			)
		}
	# now, hopefully, none of the remaining tips/nodes
		# do not have NA appearance times...	
	###########################
	# get node max ages
	if(is.null(taxaTree$node.label)){
		# no node labels...
		# so find the max age using the children on tree
		#
		nodeChildren <- lapply(prop.part(taxaTree), function(x)
			taxaTree$tip.label[x]
			)
		nodeMaxAges <- numeric()
		#
		for(i in 1:Nnode(taxaTree)){
			matchedChildren <- match(nodeChildren[[i]], appData$taxon_name)
			maxChildAges <- appData$firstapp_max_ma[matchedChildren]
			nodeMaxAges[i] <- max(maxChildAges)
			}		
	}else{
		nodeMaxAges <- appData$firstapp_max_ma[
			match(taxaTree$node.label, appData$taxon_name)
			]
		}
	# 
	####################################################
	# get four date taxon ranges for all tip taxa
	# match tip-taxa to taxa-data
		#sort based on tip labels
	tipMatches <- match(taxaTree$tip.label, appData$taxon_name)
	#
	fourDateRanges <- appData[tipMatches,
		c(firstAppTimes,lastAppTimes )]	
	# select the right tipAges based on tipTaxonDateUsed
	if(tipTaxonDateUsed == "deepestAge"){
		tipAges <- fourDateRanges$firstapp_max_ma
		}
	if(tipTaxonDateUsed == "shallowestAge"){
		tipAges <- fourDateRanges$lastapp_min_ma
		}
	###############################################
	# construct tip+node age vector
	allAges <- as.numeric(c(tipAges, nodeMaxAges))
	#print(allAges)
	#
	# check ages for NAs
	taxaNAapps <- is.na(allAges)
	#
	if(any(taxaNAapps)){
		namesAllAges <- c(taxaTree$tip.label, taxaTree$node.label)
		noappNames <- namesAllAges[taxaNAapps]
		noappNames <- paste0(noappNames, collapse = ",\n")
		stop(paste0("Some ages contain NAs: \n",
			noappNames))
		}
	##########################################
	#print(taxaTree)
	#print(allAges)
	#
	# get dated tree
	datedTree <- nodeDates2branchLengths(
		nodeDates = allAges,
		tree = taxaTree,
		allTipsModern = FALSE
		)
	############################################
	# check that the tree and its root age makes sense
	#checkRes <- checkRootTime(datedTree)			
	#################################
	if(minBranchLen>0){
		# take care of zero length branches
		datedTree <- minBranchLength(
			tree = datedTree,
			mbl = minBranchLen)
		#
		datedTree <- ladderize(datedTree)
		#
		plotName <- paste0(
			"Dated Phylogeny (",
			minBranchLen,
			" Minimum Branch Length)"
			)
	}else{
		plotName <- "Dated Phylogeny"
		}
	#
	#print("a")
	#
	############################################
	# check that the tree and its root age makes sense
	checkRootTimeRes <- checkRootTime(datedTree,
		stopIfFail = TRUE)
	#################
	if(plotTree){
		plot(datedTree, 
			main = plotName,
			show.node.label = FALSE, 
			cex = 0.5)
		axisPhylo()
		}
	###################		
	# output four date taxon ranges for all tip taxa
		# (for strap or otherwise)	
	rownames(fourDateRanges) <- taxaTree$tip.label
	datedTree$tipTaxonFourDateRanges <- fourDateRanges
	# and ranges used - return deepest and shallowest time
	datedTree$ranges.used <- fourDateRanges[,
		c("firstapp_max_ma","lastapp_min_ma")
		]
	colnames(datedTree$ranges.used) <- c("FAD", "LAD")
	#
	#print("b")
	#
	# output the original taxaDataPBDB_orig
	datedTree$taxaDataPBDB <- taxaDataPBDB_orig
	# return
	return(datedTree)
	}



processTaxonAppData <- function(dataPBDB, dropZeroOccurrenceTaxa){
	#
	lastAppTimes <- c("lastapp_min_ma","lastapp_max_ma")
	firstAppTimes <- c("firstapp_min_ma","firstapp_max_ma")
	##########################################
	# identify which taxa are extant
	isExtant <- dataPBDB$is_extant == "extant"	
	isExtinct <- dataPBDB$is_extant == "extinct"	
	hasZeroOcc <- (dataPBDB$n_occs == 0)
	isExtinctAndZeroOcc <- hasZeroOcc & isExtinct
	#
	hasNAFirstApps <- apply(
		dataPBDB[,firstAppTimes], 1,function(x) is.na(x[1]) & is.na(x[2])
		)	
	######################################################
	if(any(isExtant)){
		# replace min & max last appearance ages
			# with 0 if "is_extant" is "extant"
		dataPBDB[isExtant, lastAppTimes] <- 0
		# if both first appearance times for a taxon are NA (no fossil occurrences)
			# and its extant, put in 0 for its first appearance time
		dataPBDB[isExtant & hasNAFirstApps, firstAppTimes] <- 0
		}
	######################################################
	# 
	if(any(isExtinctAndZeroOcc) & dropZeroOccurrenceTaxa){
		#message("Some extinct nodes have zero occurrences and thus no age data - these are removed")
		dataPBDB <- dataPBDB[-which(isExtinctAndZeroOcc),]
		}
	###############################
	# test if any are still NA for all app times
	#hasNAApps <- apply(dataPBDB[,
	#		c(firstAppTimes,lastAppTimes )
	#		],
	#	1, function(x) is.na(x[1]) & is.na(x[2])
	#	)
	#dataPBDB$taxon_name[hasNAApps] 	
	#
	return(dataPBDB)
	}	

