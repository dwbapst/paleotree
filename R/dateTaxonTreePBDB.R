#' Date a Taxon-Topology from the Paleobiology Database Using Appearance Data from the API
#' 
#' The function \code{dateTaxonTreePBDB} takes a topology of  
#' The required input is a topology with tip and internal node labels corresponding to taxa in the Paleobiology Database, and a table of data (containing those same tip and node taxa) obtained from the taxa-list functionality of the Paleobiology Database's API, with appearance times output.

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

#' @return
#' Returns a dated phylogeny of class \code{phylo}, with an additional element
#' \code{$taxaDataPBDB} added containing the input \code{taxaDataPBDB}, as this might be
#' called by other functions.

#' @seealso
#' See \code{\link{getTaxaDataPBDB}}, \code{\link{makePBDBtaxonTree}},
#' and \code{\link{plotPhylopicTreePBDB}}.

#' @author David W. Bapst

#' @references
#' Peters, S. E., and M. McClennen. 2015. The Paleobiology Database
#' application programming interface. \emph{Paleobiology} 42(1):1-7.

#' @examples
#' 
#' #an example here


#' @name dateTaxonTreePBDB
#' @rdname dateTaxonTreePBDB
#' @export
dateTaxonTreePBDB <- function(
		taxaTree,
		taxaDataPBDB = taxaTree$taxaDataPBDB,
		minBranchLen = 0,
		plotTree = FALSE){
	###################################
	if(!any(colnames(taxaDataPBDB)!="lastapp_min_ma")){
		stop(paste0(
			"input taxaTree must have a taxaDataPBDB element of PBDB data",
			"generated with show=app"))
		}
	############################
	# get tip min ages
	# first replace min ages with 0 if "is_extant" is "extant"
	taxaDataPBDB$lastapp_min_ma[taxaDataPBDB$is_extant == "extant"] <- 0
	#
	# now sort based on tip labels
	tipMinAges <- taxaDataPBDB$lastapp_min_ma[
		match(taxaTree$tip.label, taxaDataPBDB$taxon_name)]
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
	if(plotTree){
		plot(datedTree, main=plotName,
			show.node.label=FALSE, cex=0.5)
		axisPhylo()
		}
	###################
	# return
	datedTree$taxaDataPBDB <- taxaDataPBDB
	return(datedTree)
	}