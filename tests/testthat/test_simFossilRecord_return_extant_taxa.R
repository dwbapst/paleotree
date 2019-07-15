test_that("simFossilRecord doesn't return extant taxa in extinct-only simulations", {

library(paleotree)

# We can set up a test to make sure that no extant taxa somehow get
# returned in many simulations with extinct-only conditioning:

set.seed(1)

res <- simFossilRecord(
    p = 0.1, 
    q = 0.1, 
    r = 0.1,
    nTotalTaxa = 10,
    nExtant = 0,
    nruns = 1000,
    plot = TRUE
    )
	
anyLive <- any(sapply(res,function(z) 
    any(sapply(z,function(x) x[[1]][5] == 1)))
    )
	
# test if any are still alive
if(anyLive){
    stop("Runs have extant taxa under conditioning for none?")
    }

expect_false(anyLive)

})
