test_that("obtainDatedPosteriorTreesMrB works", {

#library(paleotree)
#setwd("d://dave//workspace//paleotree//tests//testthat")

fileTest <- "data//retio_dating.run1.t"

MCCT <- obtainDatedPosteriorTreesMrB(
	runFile = fileTest,
	nRuns = 2, 
	burnin = 0.5,
	getFixedTimes = FALSE,
	outputTrees = "MCCT", 
	file = NULL)

expect_message(
MAP <- obtainDatedPosteriorTreesMrB(
	runFile = fileTest,
	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
	outputTrees = "MAPosteriori", 
	file = NULL)
)
		
MAPr <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
	outputTrees = "MAPriori", 
	file = NULL)

expect_message(
expect_warning(		
MaxLike <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
	outputTrees = "MaxLikelihood", 
	file = NULL)
))

# get a root age from the fixed ages for tips
	# but root age already exists so this should return an error
expect_error(
setRootAge(tree = MCCT)
)

#pull a hundred trees randomly from the posterior
hundredRandomlySelectedTrees <- obtainDatedPosteriorTreesMrB(
 	runFile = fileTest,
 	nRuns = 2, 
	burnin = 0.5, 
	getFixedTimes = TRUE,
 	getRootAges = TRUE,
	outputTrees = 100, 
	file = NULL)

testthat::skip_on_cran()
testthat::skip_on_travis()	
expect_equal_to_reference(hundredRandomlySelectedTrees) 

})