test_that("timeSliceTree tip label options work",{

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
    plot = TRUE,
    drop.extinct = FALSE
    )
# compare tip labels when we use tipLabels = "allDesc"
tree950_AD <- timeSliceTree(
    tree,
    sliceTime = 950,
    plot = TRUE,
    tipLabel = "allDesc",
    drop.extinct = FALSE
    )
    
if(all(tree950$tip.label == tree950_AD$tip.label)){
  stop("tip labels seem to be the same but should be different??")
}

})