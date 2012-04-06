\name{sampleRanges}
\alias{sampleRanges}
\title{Sampling Taxon Ranges}
\description{A function for simulating the effect of incomplete sampling of the fossil record.}
\usage{sampleRanges(taxad, r, min.taxa = 2, ranges.only = T, modern.samp = T)}
\arguments{
  \item{taxad}{A five column matrix, as output by simFossilTaxa, or a two-column matrix of per-taxon ranges}
  \item{r}{Instantaneous sampling rate per time unit}
  \item{min.taxa}{Minimum number of taxa sampled}
  \item{ranges.only}{Give time of all sampling events or ranges only?}
  \item{modern.samp}{Always perfectly sample living taxa at the present day (time=0)?}
}
\details{
Currently, this function only implements the simplest of sampling models, where sampling events occur under a homogenous rate through time, a homogenous rate across lineages and the waiting times between them are exponentially distributed across lineage*time-units (Foote, 1997). Future modifications of this function will include more sophisticated options, such as the hat model (Liow et al., 2010).

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

If min.taxa=0, the simulation may output ranges for which no taxa were ever sampled.

If modern.samp=T, then living taxa will always be sampled at least at the present day (if there are any living taxa). If not true, then it is assumed that almost nothing is known about the living taxa. This is to simulate the fact that we generally assume that is the modern day is much better sampled than the fossil record, although probably not perfectly. Presently, there is no option for increased sampling at the present-day.
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
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
taxicDivCont(rangesCont)
}
\keyword{datagen}