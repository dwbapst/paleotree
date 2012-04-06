\name{getSampRateCont}
\alias{getSampRateCont}
\title{Fit Models of Sampling Rates to Continuous-Time Taxon Ranges}
\description{Uses ML to find the best-fit parameters for models of sampling and extinction rates, given a set of continuous-time taxon ranges from the fossil record
}
\usage{
getSampRateCont(timeData, n_tbins = 1, grp1 = NA, grp2 = NA, threshold = 0.1, est_only = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeData}{Two-column matrix of per-taxon first and last occurrances in absolute continous time}
  \item{n_tbins}{Number of time bins with different sampling/extinction parameters}
  \item{grp1}{A vector of integers or characters, the same length as the number of taxa in timeData, where each taxon-wise element gives the group ID of the taxon for the respective row of timeData}
  \item{grp2}{A vector of integers or characters, the same length as the number of taxa in timeData, where each taxon-wise element gives the group ID of the taxon for the respective row of timeData}
  \item{threshold}{The smallest allowable range (measured difference in the FAD and LAD of a taxon). Ranges below this size will be treated as "one-hit" sampling events.}
  \item{est_only}{If true, function will give back a matrix of ML extinction rates and sampling probabilities per species rather than usual output (see below)}
}
\details{
This function uses maximum-likelihood solutions found by Foote (1997). These analyses are ideally applied to data from single stratigraphic section but can potentially be applicable to regional or global datasets, although the behavior of those datasets is less well understood.

getSampRateCont allows for a considerable level of versatility in terms of the degree of variation allowed among taxa in sampling rates. Essentially, this function allows taxa to be broken down into different possibly overlapping classes which have 'average' parameter values that are then combined to calcualte per-taxon parameters. For example, perhaps I think that taxa that live in a particular environment have a different characteristic sampling rate/probability, taxa of several different major clades have different characteristic sampling parameters and that there may be several temporal shifts in the characteristic extinction rate or sampling parameters. The classification IDs for the first two can be included as per-taxon vectors (of either character or integers) as grp1 and grp2 and the hypothesized number of temporal breaks can be included as the n_tbins argument. A model where taxa differ in parameters across time, clades and environments will then be fit and the AIC calculated, so that it can be compared to other models. 

By default, the simple model where all taxa belong to a single class, with a single characteristic extinction rate and a single characteristic sampling parameter, is fit to the range data.

The timebins option allows for timebins with free-floating boundaries that are not defined a priori. The boundaries between time bins with different characteristic parameters will thus be additional parameters included in the AIC calculation. If you have the prior inclination that sampling/extinction changed at a particular point in time, then seperate the taxa that originated before and after that point as two different groups and include those classifications as a grp in the arguments.

This function does not implement the finite window of observation modification for dealing with data that leads up to the recent. This is planned for a future update, however. For now, data input into this function should be for taxa that have already gone extinct by the modern and are not presently extant.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

Please check the $message element of the output to make sure that convergence occurred.
}
\value{
If est_only=T, a matrix of per-taxon sampling and extinction parameters is output.

If est_only=F (default), then the output is a list:

\item{Title}{Gives details of the analysis, such as the number and type of parameters included and the number of taxa analyzed}
\item{pars}{Maximum-likelihood parameters of the sampling model}
\item{SMax}{The maximum support (log-likelihood) value}
\item{AICc}{The second-order Akaike's Information Criterion, corrected for small sample sizes}
\item{message}{Messages output by optim(); check to make sure that model convergence occurred}

If the multi-class models are using, the element $pars will not be present, but there will be several different elements that sum the characteristic parameter components for each class. As noted in the $title, these should not be interpretated as the actual rates/probabilities of any real taxa but rather as components which must be assessed in combination with other classes to be meaningful. For example, for taxa of a given group in a given time bin, their extinction rate is the extinction rate component of that time bin times the extinction rate component of their group. Completeness estimates (Pp) will be output with these parameters as long as classes are not overlapping, as those estimates would not otherwise refer to meaningful groups of taxa.
}
\references{
Foote, M. 1997. Estimating Taxonomic Durations and Preservation Probability. Paleobiology 23(3):278-300.

Foote, M., and D. M. Raup. 1996. Fossil preservation and the stratigraphic ranges of taxa. Paleobiology 22(2):121-140.
}
\author{David W. Bapst}
\seealso{
\code{\link{getSampProbDisc}},\code{\link{sRate2sProb}},\code{\link{qsRate2Comp}},
}
\examples{
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
(SRres1<-getSampRateCont(rangesCont))
#that's all the results...
sRate<-SRres1$pars[2]
print(sRate)	#estimates that sRate=~0.4 (not too bad...)

#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins
	#AIC should be higher (less informative)
randomgroup<-sample(1:2,nrow(rangesCont),replace=TRUE)
SRres2<-getSampRateCont(rangesCont,grp1=randomgroup)
SRres3<-getSampRateCont(rangesCont,n_tbins=2)
SRres4<-getSampRateCont(rangesCont,n_tbins=3,grp1=randomgroup)
print(c(SRres1$AICc,SRres2$AICc,SRres3$AICc,SRres4$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)
}
