#' @details

#' @param taxaTree A tree with tip taxon names matching the taxon names
#' in \code{taxaData}. Probably a taxon tree estimated
#' with \code{\link{makePBDBtaxonTree}}.

#' @param taxaData A data table of taxonomic information obtained
#' using the Paleobiology Database's API for a set of taxa that
#' includes the tip taxa on \code{taxaTree}, generated with
#' parameter \code{show=app} so that appearance times are included.

#' @param minBranchLen Following dating using the appearance times taken directly
#' from the PBDB for each tip taxon and node, the tree may then be assessed with
#' the minimum branch length algorithm, as applied by \code{\link{minBranchLength}}.
#' If \code{minBranchLen = 0}, the default, this step is skipped. It may be necessary
#' to set \code{minBranchLen} higher than zero to differentiate nodes in cases
#' with poor stratigraphic congruency, so that derived taxa are the first taxa
#' observed in a group.

#' @param plot If \code{TRUE}, the resulting dated tree is plotted.

#' @return

#' @aliases

#' @seealso
#' See \code{\link{getTaxaDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhylopicTreePBDB}}.

#' @author David W. Bapst

# @references

#' @examples

#' @name dateTaxonTreePBDB
#' @rdname dateTaxonTreePBDB
#' @export
dateTaxonTreePBDB <- function(
		taxaTree,
		taxaData = taxaTree$taxaData,
		minBranchLen = 0,
		plot = FALSE){
	###################################
	taxaTree <- taxaTree
	taxaData <- taxaTree$taxaData
	if(!any(colnames(taxaData)!="lastapp_min_ma")){
		stop(paste0(
			"input taxaTree must have a taxaData element of PBDB data",
			"generated with show=app"))
		}
	############################
	# get tip min ages
	# first replace min ages with 0 if "is_extant" is "extant"
	taxaData$lastapp_min_ma[taxaData$is_extant == "extant"] <- 0
	#
	# now sort based on tip labels
	tipMinAges <- taxaData$lastapp_min_ma[
		match(taxaTree$tip.label, taxaData$taxon_name)]
	#
	# get node ages
	nodeNames<-paste0(taxaTree$node.label,collapse=",")
	apiAddressNodes <- paste0("http://paleobiodb.org/data1.2/taxa/list.txt?name=",
		nodeNames,"&show=app,parent"
		)
	# browseURL(apiAddressNodes)
	nodeData <- read.csv(apiAddressNodes ,
		stringsAsFactors = FALSE)
	# get node max ages
	nodeMaxAges <- nodeData$firstapp_max_ma[
		match(taxaTree$node.label, nodeData$taxon_name)]
	# construct tip+node age vector
	allAges <- as.numeric(c(tipMinAges,nodeMaxAges))
	##########################################
	# get dated tree
	datedTree <- nodeDates2branchLengths(
		nodeDates = allAges,
		tree = taxaTree,
		allTipsModern = FALSE)
	#
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
	#################
	if(plot){
		plot(datedTree, main=plotName,
			show.node.label=FALSE, cex=0.5)
		axisPhylo()
		}
	###################
	# return
	datedTree$taxaData <- taxaData
	return(datedTree)
	}