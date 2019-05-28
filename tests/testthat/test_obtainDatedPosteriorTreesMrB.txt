test_that("obtainDatedPosteriorTreesMrB works", {

#library(paleotree)

fileTest <- "balnibarbinae_mrB_tip-dating_08-13-18.nex.run1.t"

MCCT <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5,
	outputTrees = "MCCT", 
	file = NULL)

MAP <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
	outputTrees = "MAPosteriori", 
	file = NULL)
		
MAPr <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
	outputTrees = "MAPriori", 
	file = NULL)
		
MaxLike <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
	outputTrees = "MaxLikelihood", 
	file = NULL)

# get a root age from the fixed ages for tips
setRootAge(tree = MCCT)
setRootAge(tree = MAP)
setRootAge(tree = MAPr)
setRootAge(tree = MaxLike)

#pull a hundred trees randomly from the posterior
hundredRandomlySelectedTrees <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
 	getRootAges = TRUE,
	outputTrees = 100, 
	file = NULL)

}}