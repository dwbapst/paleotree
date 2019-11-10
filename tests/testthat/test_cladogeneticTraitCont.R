test_that("cladogeneticTraitCont functions fine", {

#library(paleotree)
set.seed(444)

record <- simFossilRecord(
	p = 0.1, q = 0.1, 
	nruns = 1,
	nTotalTaxa = c(10,50), 
	plot = TRUE)
taxa <- fossilRecord2fossilTaxa(record)
trait <- cladogeneticTraitCont(taxa)
tree <- taxa2phylo(taxa)
plotTraitgram(trait,tree,conf.int = FALSE)

#with cryptic speciation
record2 <- simFossilRecord(
	p = 0.1, q = 0.1, prop.cryptic = 0.5, 
	nruns = 1, 
	nTotalTaxa = c(10,50), 
	plot = TRUE)
taxa2 <- fossilRecord2fossilTaxa(record2)
trait2 <- cladogeneticTraitCont(taxa2)
tree2 <- taxa2phylo(taxa2)
plotTraitgram(trait2,tree2,conf.int = FALSE)

testthat::skip_on_cran()
testthat::skip_on_travis()

expect_equal_to_reference(trait, update=TRUE, 
	file = ".//references//cladogeneticTrait//trait")
expect_equal_to_reference(tree, update=TRUE, 
	file = ".//references//cladogeneticTrait//tree")

expect_equal_to_reference(trait2, update=TRUE, 
	file = ".//references//cladogeneticTrait//trait2")
expect_equal_to_reference(tree2, update=TRUE, 
	file = ".//references//cladogeneticTrait//tree2")

})
