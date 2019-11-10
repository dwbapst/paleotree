test_that("simFossilRecord is returning right output without returnAllRuns", {

# without returnAllRuns

set.seed(4)

record <- simFossilRecord(
    p = 0.1, 
    q = 0.1, 
    r = 0.1, 
    nruns = 3, 
    totalTime = c(50,100), 
    nTotalTaxa = 30, 
    nExtant = 10,
    nSamp = 15, 
    returnAllRuns = FALSE,
    print.runs = TRUE, 
    plot = TRUE
    )
	
expect_equal(length(record), 3)
}
)


test_that("simFossilRecord is returning right output for returnAllRuns", {

testthat::skip_on_cran()
testthat::skip_on_travis()
	
# with returnAllRuns

set.seed(4)

record <- simFossilRecord(
    p = 0.1, 
    q = 0.1, 
    r = 0.1, 
    nruns = 3, 
    totalTime = c(50,100), 
    nTotalTaxa = 30, 
    nExtant = 10,
    nSamp = 15, 
    returnAllRuns = TRUE,
    print.runs = TRUE, 
    plot = TRUE
    )

x <- sapply(record$rejected, divCurveFossilRecordSim)

expect_equal(length(record), 2)
expect_identical(names(record), c("accepted", "rejected"))
expect_equal(length(record$accepted), 3)

}
)