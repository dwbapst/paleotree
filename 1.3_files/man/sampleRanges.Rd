\name{sampleRanges}
\alias{sampleRanges}
\title{Sampling Taxon Ranges}
\description{A function for simulating the effect of incomplete sampling of the fossil record.}
\usage{sampleRanges(taxad, r, min.taxa = 2, ranges.only = T, modern.samp = T)}
\arguments{
  \item{taxad}{A two-column matrix of per-taxon ranges. The five-column matrix output of simFossilTaxa can also be supplied.}
  \item{r}{Instantaneous sampling rate per time unit}
  \item{min.taxa}{Minimum number of taxa sampled. The default is 2.}
  \item{ranges.only}{Give time of all sampling events or ranges only?}
  \item{modern.samp}{Probability of sampling living taxa at the present day (time=0), see below.}
}
\details{
Currently, this function only implements a simple sampling models, where sampling events occur under a rate which is homogenous through time and across lineages and the waiting times between these sampling events are exponentially distributed across lineage/time-units (Foote, 1997). Future modifications of this function will include more sophisticated options, such as the hat model (Liow et al., 2010).

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

If min.taxa=0, the simulation may output ranges for which no taxa were ever sampled.

The taxa still alive at time=0 (if there are any living taxa) will be sampled 

If modern.samp.prob=1.0 (the default), then living taxa will always be sampled at least at the present day (if there are any living taxa). If the probability is less than 1

If not true, then it is assumed that almost nothing is known about the living taxa. This is to simulate the fact that we generally assume that is the modern day is much better sampled than the fossil record, although probably not perfectly. Presently, there is no option for increased sampling at the present-day.

Conditioning on sampling some minimum number of taxa may create strange simulation results for some analyses, such as simulation analyses of birth-death processes. Set min.taxa=0 to remove this conditioning.
}
\value{
If ranges.only=T, then the output is a two-column per-taxon matrix of first and last appearances in absolute time. NAs mean the taxon was never sampled in the simulation.

If ranges.only=F, then the output is a list, where each element is a vector of sampling events the timing of sampling events, each corresponding to a different taxon in the input. Elements that are NA are unsampled taxa.
}
\references{
Foote, M. 1997. Estimating Taxonomic Durations and Preservation Probability. Paleobiology 23(3):278-300.

Liow, L. H., T. B. Quental, and C. R. Marshall. 2010. When Can Decreasing Diversification Rates Be Detected with Molecular Phylogenies and the Fossil Record? Systematic Biology 59(6):646-659.
}
\author{David W. Bapst}
\seealso{\code{\link{simFossilTaxa}},\code{\link{binTimeData}}}
\examples{
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
layout(matrix(1:2,2,))
#let's see what the 'true' diversity curve looks like in this case
taxicDivCont(taxa)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
taxicDivCont(rangesCont)
#compare the true history to what we might observe!

#let's try a more complicated model: 
}