test_that("quick example of birth-death simulation with simFossilRecord", {
    
    library(paleotree)
    set.seed(2)
    
    # quick birth-death-sampling run
        # with 1 run, 50 taxa
    
    record <- simFossilRecord(
        p = 0.1, 
        q = 0.1, 
        r = 0.1, 
        nruns = 1,
        nTotalTaxa = 50, 
        plot = TRUE
        )
  
  expect_equal_to_reference(record, 
                            update=TRUE,
                            file = ".//references//quick_simFossilRecord_record") 
  
})
