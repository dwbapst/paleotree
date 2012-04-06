\name{simPaleoTrees}
\alias{simPaleoTrees}
\title{Simulating Un-Conditioned Trees of Fossil Taxa}
\description{Creates sets of paleo-trees with minimal conditioning and with sampling, except node-times are known perfectly}
\usage{simPaleoTrees(p, q, r, ntrees=1, all.extinct=F, modern.samp=F, mintime=1, 
	maxtime=100, mintaxa=2, maxtaxa=500, drop.zlb = T, plot = F)}
\arguments{
  \item{p}{Instantaneous rate of speciation per lineage-time units}
  \item{q}{Instantaneous rate of extinction per lineage-time units}
  \item{r}{Instantaneous rate of sampling per lineage-time units}
  \item{ntrees}{Number of trees to simulate}
  \item{all.extinct}{Condition on all taxa being extinct by modern? Default is false}
  \item{modern.samp}{Should taxa extant at present be sampled perfectly at the present-day (time=0)?}
  \item{mintime}{Minimum time to run simulation}
  \item{maxtime}{Maximum time to run simulation}
  \item{mintaxa}{Minimum number of taxa allowed over entire simulation}
  \item{maxtaxa}{Maximum number of taxa allowed over entire simulation}
  \item{drop.zlb}{Should zero-length branches be dropped?}
  \item{plot}{Should data be plotted as it is simulated?}
}
\details{
Essentially a complex wrapper for simFossilTaxa, sampleRanges and taxa2phylo, this function will spit out simulated phylogenies of fossil taxa where the divergence times are perfectly known, only sampled taxa are included on the tree and tips are located at the last observed time for the species (the apparent time of extinction, except for living taxa).

simPaleoTrees essentially uses simFossilTaxa with no.cond=T and other minimal conditioning, to get as unbiased a sample of simulations as possible (without exceeding the maximum constraints). This is useful for birth-death analyses, although a number of options available in simFossilTaxa are thus unavailable in simPaleoTrees. By default, there is no conditioning on the number of extant taxa, living taxa are sampled perfectly at time 0 and zero-length branches are dropped. Unlike simPaleoTrees, you cannot condition on a certain number of extant taxa, only whether they are allowed or not (via all.extinct). There are also no options relating to speciation modes: taxa are only simulated under budding cladogenesis.

Because the divergence times are known perfectly, yet tips are at the apparent time of extinction and unsampled taxa are dropped, one should not use the output of this analysis except for very specialized simulation analyses. The results are probably not anything like real data, at least in most qualities.
}
\value{
Output is an object of class multiphylo containing the simulated phylogenies, unless ntrees=1 in which case the output is a phylogeny of class 'phylo'. Additionally, each of these simulated phylogenies will have the original simulated taxa data (from simFossilTaxa) and sampled ranges (from sampleRanges) attached as the elements $taxa and $ranges.
}
\author{David W. Bapst}
\seealso{
\code{\link{simFossilTaxa}},\code{\link{taxa2phylo}},\code{\link{sampleRanges}},
}
\examples{
dev.new()
set.seed(444)
#simulate trees conditioned to have no living descendants
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=10,all.extinct=TRUE,maxtime=100,plot=TRUE)
#number of tips
sapply(trees,Ntip)

#simulate trees conditioned to (possibly) have living descendants and perfect sampling at modern
trees<-simPaleoTrees(p=0.1,q=0.1,r=0.5,ntrees=10,all.extinct=FALSE,maxtime=100,modern.samp=TRUE,plot=TRUE)
#number of tips
sapply(trees,Ntip)
}
