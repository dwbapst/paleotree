

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
record <- simFossilRecord(
	p = 0.1, q = 0.1, prop.cryptic = 0.5, 
	nruns = 1, 
	nTotalTaxa = c(10,50), 
	plot = TRUE)
taxa <- fossilRecord2fossilTaxa(record)
trait <- cladogeneticTraitCont(taxa)
tree <- taxa2phylo(taxa)
plotTraitgram(trait,tree,conf.int = FALSE)
