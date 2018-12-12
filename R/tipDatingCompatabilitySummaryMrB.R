#' Get the Compatibility Summary Topology From a Tip-Dating Analysis with MrBayes
#' 
#' This function is designed to avoid methodological
#' issues with getting sensible consensus summary topologies from posteriors samples of
#' tip-dated, sampled-ancestor trees output by Mr Bayes. This function will obtain samples of posterior trees,
#' from external files, remove the specified burn-in,
#' and output an undated summary tree of clades (splits) indicated on the output tree, as a particular posterior 
#' probability threshold. Posterior probabilities may be appended to the nodes of the output phylogeny.
#' This function should be used for examining topological variation in the posterior.

#' @details
#' This function is most useful for dealing with dating analyses in MrBayes, particularly when tip-dating
#' a tree with fossil taxa, as the half-compatibility and all-compatibility summary trees offered by the
#' '\code{sumt}' command in MrBayes can have issues properly portraying summary trees from such datasets.
#' 
#' Summary topologies calculated with some tip-dating software environments, such as MrBayes,
#' can be subject to strange and uninterpretable methodological artifacts
#' as the methods use attempt to present summary topologies with
#' branch lengths. Many of these algorithms as currently implemented
#' cannot handle the two-degree nodes or zero-length branches
#' that arise from having sampled ancestors. Users looking to summarize a tip-dating analysis cannot
#' easily calculate a dated summary: if they want a dated tree, they \emph{must} examine a
#' single tree from the posterior (either randomly selected or chosen based on some criteria
#' such as marginal likelihood, posterior probability, etc). However, if our main interest is
#' the unscaled evolutionary closeness of taxonomic units without reference to time,
#' then it is sufficient to examine a summary of the topological variation over our posterior.
#' 

#' @inheritParams obtainDatedPosteriorTreesMrB

#' @param compatibilityThreshold The posterior probability threshold (between 1 and zero, post-burn-in)
#' that a node must satisfy to appear on the output summary tree. The default is 0.5, making the
#' trees output half-compatibility trees (summary topologies), similar to the majority-rule
#' consensus commonly used in maximum parsimony analyses.  The value cannot be lower than 0.5 due to
#' current technical constraints, and the need for an R function that iteratively ranks possible splits to
#' be included in a consensus, as the consensus is calculated. Currently, if a clade frequency threshold
#' given (argument \code{p}) to  \code{ape} function \code{\link[ape]{consensus}}, which is used
#' internally by \code{halfCompatTree}, all nodes above that compatibility threshold, even splits
#' which are contradictory, will be included on the output tree, often resulting in uninterpretable output.

#' @note
#' Consensus trees that combine clades found different trees in the same tree sample
#' may inadvertently combine clades that are not found on any of the actual trees
#' sampled in the posterior, and may be quite far from the posterior trees
#' as sampled in multivariate tree-space. This is a standard criticism leveled at consensus-type summary trees,
#' except for the strict consensus (equivalent here to if a user tried \code{compatibilityThreshold = 1}). However,
#' post-burnin posterior tree samples often sample (and thus contain) a considerable range of
#' tree-space within them, and thus the strict consensus (a total compatibility tree?)

#' @return
#' A single, undated summary tree, containing those clades (splits) found in greater frequency in
#' the post-burnin posterior tree sample more than the value of \code{compatibilityThreshold}, of
#' class \code{phylo}. If \code{labelPostProb = TRUE}, nodes will be labeled with the posterior probability of
#' the respective clade.

#' @seealso
#' See function \code{\link{obtainDatedPosteriorTreesMrB}} for additional
#' ways of processing and evaluating trees from MrBayes posterior samples.
#' 
#' Summary trees are estimated using the function \code{\link[ape]{consensus}}
#' in package \code{ape}.

#' @author David W. Bapst

#' @examples
#' \dontrun{
#' #pull post-burnin trees from the posterior
#'       # and get the half-compatibility summary (majority-rule consensus)
#'       # by setting 'compatibilityThreshold = 0.5'
#' 
#' halfCompatTree <- tipDatingCompatabilitySummaryMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, burnin = 0.5, 
#'  	compatibilityThreshold = 0.5,
#'  	labelPostProb = TRUE
#'  	)
#' 
#' # let's try plotting it with posterior probabilities as node labels
#' plot(halfCompatTree)
#' nodelabels(halfCompatTree$node.label)
#'  
#' }
#' 

#' @name tipDatingCompatabilitySummaryMrB
#' @rdname tipDatingCompatabilitySummaryMrB
#' @export
tipDatingCompatabilitySummaryMrB <- function(
		runFile,
		nRuns = 2,
		burnin = 0.5,
		compatibilityThreshold = 0.5,
		labelPostProb = TRUE
		){
	########################################################
	# CHECKS
	if(compatibilityThreshold<0.5){
		stop("compatibilityThreshold < 0.5: Currently cannot calculate compatibility trees containing nodes found on less than half of all trees")
		}
	## CAN WE DO THIS?
	#If compatibilityThreshold = 0 is used, 
		# the tree is effectively a total-compatibility tree,
		# fully bifurcating, with relationships resolved in order
		# of posterior probabilities (post-burnin).
	# 'extended majority consensus tree' or 'all-compatible-nodes (allcompat) tree'
	#
	##########################################################
	# #pull post-burnin trees from the posterior
	postBurnInTrees <- obtainDatedPosteriorTreesMrB(
		runFile = runFile,
		nRuns = nRuns, 
		burnin = burnin, 
		getFixedTimes = FALSE,
		getRootAges = FALSE,
		originalNexusFile = NULL,
		outputTrees = "all",
		file = NULL)
	# remove branch lengths
	for(i in 1:length(postBurnInTrees)){
		postBurnInTrees[[i]]$edge.length<-NULL
		}
	#
	# get compatibility tree for a given p using ape function consensus
	finalTree <- consensus(postBurnInTrees, p = compatibilityThreshold)
	#
	# get the posterior probabilities
	if(labelPostProb){ 
		# label tree with posterior probabilities as node labels
		finalTree <- getPosteriorProbabiities(finalTree,postBurnInTrees)
		}
	# return final compatibility tree
	return(finalTree)
	}
