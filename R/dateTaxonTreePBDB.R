#' @details

#' @inheritParams

#' @param

#' @return

#' @aliases

#' @seealso

#' @author David W. Bapst

#' @references

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
			"generated using"))
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