\name{taxa2phylo}
\alias{taxa2phylo}
\title{Convert taxon data into Phylogeny}
\description{Convert temporal and ancestor-descendant relationships of taxa into a time-scaled phylogeny}
\usage{taxa2phylo(taxad, obs_time = NULL, plot = F)}
\arguments{
  \item{taxad}{A five-column matrix of taxonomic data, as output by simFossilTaxa}
  \item{obs_time}{Per-taxon times of observation; if NULL, the LADs (column 4) in taxad2 are used}
  \item{plot}{Plot the resulting phylogeny?}
}
\details{
As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

All data relating to when static morpho-taxa appear or disappear in the record is lost; branching points will be the actual time of speciation, which (under budding) will often be in the middle of the temporal range of a taxon.
}
\note{
Please, DO NOT use this function to time-scale a real tree for a real dataset. It assumes you know the divergence/speciation times of the branching nodes perfectly, which is almost impossible given the nature of the fossil record. Use timePaleoPhy or srcTimePaleoPhy instead.

DO use this function when doing simulations and you want to make a tree of the 'true' history, such as for simulating trait evolution along phylogenetic branches.
}
\value{
The resulting phylogeny with branch lengths is output as an object of class phylo. This function will output trees with the element $root.time, which is the time of the root divergence in absolute time.
}
\author{David W. Bapst}
\seealso{
\code{\link{simFossilTaxa}},\code{\link{taxa2cladogram}},
}
\examples{
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
tree<-taxa2phylo(taxa)
phyloDiv(tree)
#now a phylogeny with tips placed at the apparent time of extinction
rangesCont<-sampleRanges(taxa,r=0.5)
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
phyloDiv(tree,drop.ZLB=FALSE)
}