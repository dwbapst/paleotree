
test_that("test involving examing the retained rejected simulations", {

    testthat::skip_on_cran()
  
    # Retaining Rejected Simulations
    
    # sometimes we might want to look at all the simulations
        # that don't meet acceptability criteria
    
    # In particular, look at simulated clades that go extinct
        # rather than surviving long enough to satisfy 
        # conditioning on temporal duration.
    
    # Let's look for 10 simulations with following conditioning:
        # that are exactly 10 time-units in duration
        # that have between 10 and 30 total taxa
        # and have 1 to 30 extant taxa after 10 time-units
    
    set.seed(4)
    
    record <- simFossilRecord(
        p = 0.1, 
        q = 0.1, 
        r = 0.1, 
        nruns = 10, 
        totalTime = 10, 
        nTotalTaxa = c(10,30), 
        nExtant = c(1,30),
        returnAllRuns = TRUE,
        print.runs = TRUE, 
        plot = TRUE
        )
    
    # when returnAllRuns = TRUE, the length of record is 2
        # named 'accepted' and 'rejected'
    
    # all the accepted runs (all 10) are in 'accepted'
    expect_equal(length(record$accepted), expected = 10)
    
    # how many taxa are in each rejected simulation run?
    totalTaxa_rej <- sapply(record$rejected, length)
    
    # plot as a histogram
    hist(totalTaxa_rej)
    # a very nice exponential distribution...
    
    # plot the rejected simulation with the most taxa
    
    divCurveFossilRecordSim(
         fossilRecord = record$rejected[[
               which(max(totalTaxa_rej) == totalTaxa_rej)[1]
               ]]
         )
    
    # we can plot all of these too...
    result <- sapply(record$rejected, 
         divCurveFossilRecordSim)
    
    # let's look at the temporal duration of rejected clades
    
    # need to write a function
    getDuration <- function(record){
         taxa <- fossilRecord2fossilTaxa(record)
         maxAge <- max(taxa[,"orig.time"], na.rm = TRUE)
         minAge <- min(taxa[,"ext.time"], na.rm = TRUE)
         cladeDuration <- maxAge - minAge
         return(cladeDuration)
         }
    
    # all the accepted simulations should have
           # identical durations (10 time-units)
    expect_true(
        all(sapply(record$accepted, getDuration) == 10)
        )
    
})