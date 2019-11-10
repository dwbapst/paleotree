test_that("timeSliceTree shifts root age when dropping extinct taxa",{

library(paleotree)

set.seed(444)
record <- simFossilRecord(
    p = 0.1, q = 0.1, nruns = 1,
    nTotalTaxa = c(30,40), 
    nExtant = 0)
taxa <- fossilRecord2fossilTaxa(record)
tree <- taxa2phylo(taxa)

tree950 <- timeSliceTree(
    tree,
    sliceTime = 950,
    drop.extinct = FALSE
    )

# compare tip labels when we use tipLabels = "allDesc"
tree950_NoExtinct <- timeSliceTree(
    tree,
    sliceTime = 950,
    drop.extinct = TRUE
    )
    
if(tree950$root.time == tree950_NoExtinct$root.time){
  stop("$root.age did not shift when extinct taxa were dropped from the time-sliced tree")
}

})