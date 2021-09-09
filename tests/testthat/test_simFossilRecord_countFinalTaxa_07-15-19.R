test_that("simFossilRecord counts taxa in final output correctly", {

library(paleotree)

set.seed(4444)

res <- simFossilRecord(
	p = 0.1, 
	q = 0.1, 
	r = 0.1,
	nTotalTaxa = 10,
	nExtant = 0,
	nruns = 1,
	plot = TRUE
	)

testthat::skip_on_cran()
testthat::skip_on_ci()

expect_equal_to_reference(res, update =TRUE,
	file = ".//references//countFinalTaxa//res")
	
})