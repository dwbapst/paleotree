\name{sampleRanges}
\alias{sampleRanges}
\title{Sampling Taxon Ranges}
\description{A function for simulating the effect of incomplete sampling of the fossil record.}
\usage{sampleRanges(taxad, r, alpha = 1, beta = 1, rTimeRatio = 1, modern.samp.prob = 1, 
	min.taxa = 2, ranges.only = TRUE, minInt = 0.01, alt.method = FALSE, plot = TRUE)}
\arguments{
  \item{taxad}{A two-column matrix of per-taxon ranges. The five-column matrix output of simFossilTaxa can also be supplied.}
  \item{r}{Instantaneous average sampling rate per lineage time units; given as a vector of length one or length equal to the number of taxa}
  \item{alpha}{Alpha parameter of beta distribution; given as a vector of length one or length equal to the number of taxa}
  \item{beta}{Beta parameter of beta distribution; given as a vector of length one or length equal to the number of taxa}
  \item{rTimeRatio}{Ratio of most recent sampling rate over earliest sampling rate; given as a vector of length one or length equal to the number of taxa}
  \item{modern.samp.prob}{Probability of sampling living taxa at the present day (time=0), see below.}
  \item{min.taxa}{Minimum number of taxa sampled. The default is 2.}
  \item{ranges.only}{Give time of all sampling events as a list or taxon ranges only?}
  \item{minInt}{Minimum interval size used for simulating complex models}
  \item{alt.method}{Use the alternative method of discretizing time even if a simple model of sampling is being simulated?}
  \item{plot}{If true, plots the sampling models for each taxon against time}
}
\details{
This function implements a range of sampling models in continuous time. Be default, sampling is simulated under the simplest model, where each sampling event occurs as a Poisson process under a instantaneous sampling rate (r) which is homogenous through time and across lineages (Foote, 1997). Under this model, the waiting times to sampling events are exponentially distributed, with an average waiting time of 1/r. This useful property allows sampling to be rapidly simulated for many taxa under this simple model in sampleRanges, by repeatedly drawing waiting times between sampling events. This is the model that is run as long as alpha, beta and rTimeRatio are set to 1.

In addition to this simple model, sampleRanges also can consider a range of additional models, in particular including the hatP and incP options of Liow et al. (2010). To describe the behavior of these models, users alter the default values for alpha, beta and rTimeRatio. These parameters, and r, can either be a single value which describes the behavior of the entire dataset or vector, of same length as the number of taxa, which describes the per-taxon value.  When any rTimeRatio, alpha or beta value is not equal to one, then the sampling rate will vary across duration of a taxon's temporal range. In general, setting alpha and beta equal to a value above 2 produce "hat" or bell-shaped curves where sampling rates peak at the midpoint of taxon ranges, while setting them unequal will produce assymmetric bell curves as in the beta function (Liow et al., 2010; Liow et al. set alpha=beta=4). rTimeRatio is the ratio of the sampling rate of the latest/most recent time divided by the earliest/oldest time.

The input r values will be interpretated differently based on whether one r value or per-taxon values were used. If one value was input, then it is assumed that r represent the grand mean r for entire dataset for purporses of time-varying r, such that if rTimeRatio is not equal to 1, taxa near the end and start of the dataset will have very different per-taxon mean sampling rate. If per-taxon values of r were input, then each r is consider the per-taxon mean sampling rate. These will not be changed, but any within-lineage variation is distributed so that the mean is still the input per-taxon value. This also changes the interpretation of rTimeRatio, such that when a single r value and rTimeRatio is given, it is assumed the ratio describes the change in sampling rates from the start of the dataset to the end, while if multiple values are given for either r or rTimeRatio will instead see the value as describing the ratio at the first and last times of each taxon. 

These distinctions between models are important: it means that all models simulated in sampleRanges are effectively nested inside the most general model which contains alpha, beta and rTimeRatio.

Because sampling rate varies over time within lineages under any of these more complex models, sampling events cannot be quickly simulated as waiting times pulled from an exponential model. Instead, the taxon durations are discretized into a large number of small time intervals of length minInt (see above; minInt should be small enough that only one event could feasibly happen per interval) and the probability of an event occurring within each interval calculated and used to simulate sampling events. In general, this method is slower but otherwise comparable to the quicker waiting times method. See examples for a small test of this.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

If min.taxa=0, the simulation may output ranges for which no taxa were ever sampled.

If modern.samp.prob=1.0 (the default), then living taxa will always be sampled at least at the present day (if there are any living taxa). If the probability is less than 1, they will be sampled with that probability at the modern day.

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
layout(1)
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

#let's try more complicated models!

#a pull-to-the-recent model with x5 increase over time similar to Liow et al.'s  incP
layout(1:2)
rangesCont1<-sampleRanges(taxa,r=0.5,rTimeRatio=5,plot=T)
taxicDivCont(rangesCont1)

#a hat-shaped model
layout(1:2)
rangesCont1<-sampleRanges(taxa,r=0.5,alpha=4,beta=4,plot=T)
taxicDivCont(rangesCont1)

#a combination of these
layout(1:2)
rangesCont1<-sampleRanges(taxa,r=0.5,alpha=4,beta=4,rTimeRatio=5,plot=T)
taxicDivCont(rangesCont1)

\dontrun{
#simulate a model where sampling rate evolves under brownian motion
tree<-taxa2phylo(taxa,obs=taxa[,3])
sampRateBM<-rTraitCont(tree)
sampRateBM<-sampRateBM-min(sampRateBM)
layout(1:2)
rangesCont1<-sampleRanges(taxa,r=sampRateBM,plot=T)
taxicDivCont(rangesCont1)

#evolving sampling rate, hat model and pull of the recent
layout(1:2)
rangesCont1<-sampleRanges(taxa,r=sampRateBM,alpha=4,beta=4,rTimeRatio=5,plot=T)
taxicDivCont(rangesCont1)

#the simpler model is simulated by pulling waiting times from an exponential
#more complicated models are simulated by discretizing time into small intervals
#are these two methods comparable?
#let's look at the number of taxa sampled under both methods
summary(replicate(100,sum(!is.na(sampleRanges(taxa,r=0.5,alt.method=F)[,1]))))
summary(replicate(100,sum(!is.na(sampleRanges(taxa,r=0.5,alt.method=T)[,1]))))
#they look pretty similar!
}
}