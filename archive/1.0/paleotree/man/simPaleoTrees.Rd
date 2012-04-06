\name{simPaleoTrees}
\alias{simPaleoTrees}
\title{Simulating Un-Conditioned Trees of Fossil Taxa}
\description{Creates sets of paleo-trees with minimal conditioning and with sampling, except node-times are known perfectly}
\usage{simPaleoTrees(p, q, r, ntrees, nExtant = NULL, maxtime=1000, maxtaxa=3000, modern.samp = T, drop.zlb = T, plot = F)}
\arguments{
  \item{p}{Instantaneous rate of speciation per lineage-time units}
  \item{q}{Instantaneous rate of extinction per lineage-time units}
  \item{r}{Instantaneous rate of sampling per lineage-time units}
  \item{ntrees}{Number of trees to simulate}
  \item{nExtant}{Maximum number of extant lineages}
  \item{maxtime}{Maximum time to run simulation}
  \item{maxtaxa}{Maximum number of taxa allowed over entire simulation}
  \item{modern.samp}{Should taxa extant at present be sampled perfectly at the present-day (time=0)?}
  \item{drop.zlb}{Should zero-length branches be dropped?}
  \item{plot}{Should data be plotted as it is simulated?}
}
\details{
Essentially a complex wrapper for simFossilTaxa, sampleRanges and taxa2phylo, this function will spit out simulated phylogenies of fossil taxa where the divergence times are perfectly known, only sampled taxa are included on the tree and tips are located at the last observed time for the species (the apparent time of extinction, except for living taxa).

The defaults are set such that there is as minimal conditioning as possible, to get as unbiased a sample of simulations as possible. There is no conditioning on the number of extant taxa, living taxa are sampled perfectly at time 0 and zero-length branches are dropped. maxtaxa and maxtime are as high as possible without decreasing performance.

Because the divergence times are known perfectly, yet tips are at the apparent time of extinction and unsampled taxa are dropped, one should not use the output of this analysis except for very specialized simulation analyses. The results are probably not anything like real data at least in most qualities.
}
\value{
Output is an object of class multiphylo containing the simulated phylogenies.
}
\author{David W. Bapst}
\seealso{
\code{\link{simFossilTaxa}},\code{\link{taxa2phylo}},\code{\link{sampleRanges}},
}
\examples{
#simulate trees with no living descendants
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=2,nExtant=0,modern.samp=TRUE,drop.zlb=TRUE,plot=TRUE)
sapply(trees,Ntip)
#simulate trees with living descendants
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=2,nExtant=10,modern.samp=TRUE,drop.zlb=TRUE,plot=TRUE)
sapply(trees,Ntip)
}
\keyword{datagen}
