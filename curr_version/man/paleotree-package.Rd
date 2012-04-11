\name{paleotree-package}
\alias{paleotree-package}
\alias{paleotree}
\docType{package}
\title{paleotree: Paleontological and Phylogenetic Analyses of Evolution}
\description{
Analyzes, time-scales and simulates phylogenies of extinct/fossil lineages. Also plots diversity curves for stratigraphic range data and phylogenies, including combinations of these two data types.
}
\details{
\tabular{ll}{
Package: \tab paleotree\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2012-02-23\cr
License: \tab  GPL (>= 2)\cr
}
This package contains functions for analyzing sampling rates given ranges of fossil taxa, in both continuous and discrete time, functions for time-scaling phylogenies of fossil taxa and functions for simulating the fossil record in both a taxic and phylogenetic varieties and converting between these.
}
\author{
David W. Bapst

Maintainer: David W. Bapst <dwbapst@gmail.com>
}
\references{
Bapst, in prep. paleotree: an R package for paleontological and phylogenetic analyses of evolution. To be submitted to Methods in Ecology and Evolution.

Bapst, in prep. Time-scaling Trees of Fossil Taxa. To be submitted to Paleobiology.
}
\seealso{
This package relies extensively on the phylogenetic toolkit offered by the \code{\link[ape:ape-package]{ape}} package.
}
\examples{
##Simulate some fossil ranges with simFossilTaxa()
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa)

#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
layout(matrix(1:2,,2))
taxicDivCont(rangesCont)

#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)
#compare to the continuous time diversity curve

#taxa2phylo assumes we know speciation events perfectly... what if we don't?
#first, let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy() using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#plot diversity curve 
phyloDiv(ttree,drop.ZLB=TRUE)

dev.new()
#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#plot diversity curve 
dev.new()
phyloDiv(ttree)
}
