#' Place a Non-Ultrametric Tree of Fossil Taxa on Absolute Time
#' 
#' This function uses a table of fixed dates for operational-taxon-units (tip taxa) to calculate the absolute
#' age of the root divergence for a tree with branch lengths, and then appends this root age to the tree
#' as a \code{$root.time} element, and then outputs the tree. Function \code{setRootAges} is a wrapper for
#' \code{setRootAge} for use with multiple trees in a "multiPhylo" object. This function was mainly written for
#' dealing with trees of extinct taxa dated in units of absolute time from Bayesian analyses, such as with MrBayes, 
#' with trees scaled to time units by functions such as \code{\link{obtainDatedPosteriorTreesMrB}}.

#' @details
#' Trees of fossil taxa come with one issue rarely encountered by those dealing with molecular
#' phylogenies: the absolute timing of when tips and divergences is not certain. With the vast
#' majority of molecular phylogenies, it can be assumed the youngest tips occur at 0 time -- i.e.,
#' the modern. This knowledge gives the tree an anchor for figuring out the absolute timing of events.
#' Many programs and other software designed for depicting and analyzing phylogenetic hypotheses
#' assumes such an apparent absolute time-scale (in R and elsewhere). A tree of Paleozoic brachiopods that include no
#' extant members has no such anchor at time = 0, and such a default assumption in available
#' software can be misleading. The $root.time protocol is intended to grant this
#' absolute time-scale to a dated tree of fossil taxa, and is appended by most of the
#' dating functions in package paleotree. However, trees dated by other approaches, such as via tip-dating in
#' programs such as MrBayes and BEAST, will not have $root.time elements when read into R.
#' 

#' @param tree A phylogeny with branch lengths of type "phylo".

#' @param trees A list of type "multiPhylo" consisting of multiple phylogenetic trees with branch lengths.

#' @param fixedAges A table of fixed ages for tip taxa, generally as a dataframe where the
#' first column is of type character, and the second column is of type numeric.  = Such a table is automatically
#' generated as an attribute of the output from \code{\link{obtainDatedPosteriorTreesMrB}}, 
#' when argument \code{getFixedTimes = TRUE}.

#' @return
#' The input tree is output, with a new \code{$root.time} element.

#' @seealso
#' \code{setRootAges} is designed to work by default with trees on relative
#' time-scales dated by \code{\link{obtainDatedPosteriorTreesMrB}}, particularly
#' when the argument with \code{getFixedTimes = TRUE} is used to obtain fixed tip
#' ages for anchoring the tree against an absolute time-scale.

#' @author David W. Bapst

#' @examples
#'
#' set.seed(444)
#' tree <- rtree(10)
#' tipAges <- cbind(c("t1","t2"), c(15,10))
#' 
#' absTimeTree <- setRootAge(tree = tree,tipAges)
#' 
#' plot(absTimeTree)
#' axisPhylo()
#' 


#' @name setRootAge
#' @rdname setRootAge
#' @export
setRootAge <- function(tree,fixedAges = NULL){
	# function for scaling posterior trees to absolute time (get root.time)
	#
	# For all trees to be comparable, we will use the $root.time convention from paleotree (Bapst, 2012)
	if(!is(tree,"phylo")){
		stop("tree must be of type 'phylo'")
		}
	if(is.null(fixedAges)){
		if(is.null(attr(tree,"fixedTable"))){
			stop("fixedAges must be supplied")
		}else{
			fixedAges <- attr(tree,"fixedTable")
			}
		}
	if(!is.null(tree$root.time)){
		stop("why does tree already have a $root.time element?? Remove to run this function")
		}
	fixedTaxa <- as.character(fixedAges[,1])
	fixedAges <- as.numeric(fixedAges[,2])
	taxaTree <- tree$tip.label
	# drop unshared taxa
	missingAge <- sapply(fixedTaxa,function(x) all(x != taxaTree))
	if(sum(missingAge) == length(fixedTaxa)){
		stop("None of the taxa in fixedAges found as OTU tip labels on tree")
		}	
	fixedAges <- fixedAges[!missingAge]
	fixedTaxa <- fixedTaxa[!missingAge]
	# get first taxon at youngest age
	youngest <- which(fixedAges == min(fixedAges))[1]
	youngDate <- fixedAges[youngest]
	youngTipDepth <- node.depth.edgelength(tree)[1:Ntip(tree)][taxaTree == fixedTaxa[youngest]]
	tree$root.time <- youngTipDepth+youngDate	
	if(length(tree$root.time)>1){
		stop("Multiple root.time elements calculated??!")}
	return(tree)
	}

#' @rdname setRootAge
#' @export
setRootAges <- function(trees,fixedAges = NULL){
	if(!is(trees,"multiPhylo")){
		stop("trees must be of type multiPhylo")
		}
	if(!is.list(trees)){
		stop("trees must be of type multiPhylo")
		}
	#
	for(i in 1:length(trees)){
		tree <- trees[[i]]
		if(!is.null(attr(trees,"fixedTable"))){
			attr(tree,"fixedTable") <- attr(trees,"fixedTable")
			}
		trees[[i]] <- setRootAge(tree = tree,fixedAges = fixedAges)
		}
	#
	return(trees)
	}
