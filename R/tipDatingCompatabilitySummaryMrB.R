#' Get the Compatibility Summary Topology From a Tip-Dating Analysis with MrBayes
#' 
#' The nature of most tip-dating analyses with MrBayes and BEAST2 will be largely incompatible with many of
#' the methods for calculating summary topologies, as these attempt to present summary topologies with
#' branch lengths -- and many of these algorithms cannot handle the two-degree nodes or zero-length branches
#' that arise from having sampled ancestors. Users looking to summarize a tip-dating analysis cannot
#' easily calculate a dated summary: if they want a dated tree, they must examine 
#' single tree from the posterior, either randomly selected or chosen based on some criteria
#' such as marginal likelihood, posterior probability, etc. However, if our main interest is examining
#' summaries of the variation in the topology over our posterior, we can go much further.



# tipDatingCompatabilitySummaryMrB

#' @inheritParams obtainDatedPosteriorTreesMrB


#' @example
#' /dontrun{
#' #pull post-burnin trees from the posterior
#' halfCompatTree <- tipDatingCompatabilitySummaryMrB(
#'  	runFile = "C:\\myTipDatingAnalysis\\MrB_run_fossil_05-10-17.nex.run1.t",
#'  	nRuns = 2, burnin = 0.5, 
#'  	compatibilityThreshold = 0.5,
#'  	labelPostProb = TRUE
#'  	)
#' }
#' 

""

#' @param compatibilityThreshold The posterior probability threshold (between 1 and zero, post-burn-in)
#' that a node must satisfy to appear on the output summary tree. The default is 0.5, making the
#' trees output half-compatibility trees (summary topologies), similar to a majority-rule
#' consensus in maximum parsimony analyses. 

## CAN WE DO THIS?
#If 0 is used, the tree is effectively a total-compatibility tree, fully bifurcating, with relationships resolved in order
# of poterior probabilities (post-burnin), such that the 


tipDatingCompatabilitySummaryMrB <- function(
	runFile,
	nRuns = 2,
	burnin = 0.5,
	compatibilityThreshold = 0.5,
	labelPostProb = TRUE
	){

	
	

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
	
	# get the 
	
	# get the posterior probabilities
	if(labelPostProb){ 
	
	getPosteriorProbabiities <- function(tree,postBurninTrees)
	
	# label tree with posterior probabilities as node labels
	}
	# return final compatibility tree
	return()
	