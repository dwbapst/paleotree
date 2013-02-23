\name{paleotree-package}
\alias{paleotree-package}
\alias{paleotree}
\docType{package}
\title{paleotree: Paleontological and Phylogenetic Analyses of Evolution}
\description{
Analyzes, time-scales and simulates phylogenies of extinct/fossil lineages, along with calculation of diversity curves. Also fits likelihood models to estimate sampling rates from stratigraphic ranges.
}
\details{
\tabular{ll}{
Package: \tab paleotree\cr
Type: \tab Package\cr
Version: \tab 1.8\cr
Date: \tab 02-22-13\cr
License: \tab  GPL (>= 2)\cr
}
This package contains functions for analyzing sampling rates given ranges of fossil taxa, in both continuous and discrete time, functions for time-scaling phylogenies of fossil taxa and functions for simulating the fossil record in both taxic and phylogenetic varieties.
}
\author{
David W. Bapst

Maintainer: David W. Bapst <dwbapst@gmail.com>
}
\references{
Bapst, D.W. 2012. paleotree: an R package for paleontological and phylogenetic analyses of evolution. \emph{Methods in Ecology and Evolution}. 3: 803-807. doi: 10.1111/j.2041-210X.2012.00223.x

Bapst, in prep. Time-scaling Trees of Fossil Taxa. To be submitted to \emph{Paleobiology}.
}
\seealso{
This package relies extensively on the phylogenetic toolkit and standards offered by the \code{\link[ape:ape-package]{ape}} package.
}
\examples{
##Simulate some fossil ranges with simFossilTaxa
set.seed(444);
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa)

#simulate a fossil record with imperfect sampling with sampleRanges
rangesCont <- sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
layout(1:2)
taxicDivCont(rangesCont)
#Now let's use binTimeData to bin in intervals of 10 time units
rangesDisc <- binTimeData(rangesCont,int.length=10)
#plot with taxicDivDisc
taxicDivDisc(rangesDisc)
#compare to the continuous time diversity curve

#taxa2phylo assumes we know speciation events perfectly... what if we don't?
#first, let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram <- taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy using the continuous range data
ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#plot diversity curve
phyloDiv(ttree,drop.ZLB=TRUE)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)
}
