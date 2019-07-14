test_that("simFossilRecord doesn't return extinct taxa in extant-only simulations", {

library(paleotree)

# We can set up a test to make sure that no extant taxa somehow get
# returned in many simulations with extinct-only conditioning:

set.seed(1)

res <- simFossilRecord(
    p = 0.1, 
    q = 0.1, 
    r = 0,
    nTotalTaxa = 10,
    nExtant = 10,
    nruns = 1000,
    plot = TRUE
    )
	
anyDead <- any(sapply(res,function(z) 
    any(sapply(z,function(x) x[[1]][5] != 1)))
    )
	
# test if any are still alive
if(anyLive){
    stop("Runs have extinct taxa under conditioning for extant only?")
    }

})
