\name{getSampProbDisc}
\alias{getSampProbDisc}
\title{Fit Models of Sampling Probability to Discrete-Interval Taxon Ranges}
\description{Uses ML to find the best-fit parameters for models of sampling probability and extinction rate, given a set of discrete-interval taxon ranges from the fossil record. This function can fit models where there are different groupings of taxa with different parameters and different free-floating time intervals with different parameters.}
\usage{
getSampProbDisc(timeData, n_tbins = 1, grp1 = NA, grp2 = NA, est_only = F, throttle = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeData}{A 2 column matrix with the first and last occurances of taxa given in relative time intervals OR the output of binTimeData}
  \item{n_tbins}{Number of time bins with different sampling/extinction parameters}
  \item{grp1}{A vector, the same length as the number of taxa in timeData, each element is a different identified for the group ID of included taxa}
  \item{grp2}{A vector, the same length as the number of taxa in timeData, each element is a different identified for the group ID of included taxa}
  \item{est_only}{If true, function will give back a matrix of ML extinction rates and sampling probabilities per species rather than usual output (see below)}
  \item{throttle}{Controls the tradeoff between run-time and accuracy}
}
\details{
This function uses maximum-likelihood solutions found by Foote (1997). These analyses are ideally applied to data from single stratigraphic section but can potentially be applicable to regional or global datasets (Foote and Raup, 1996, tested the method using Alroy's North American mammal data), although their behavior for those datasets is less well understood.

This function allows for a considerable level of versatility in terms of the degree of variation allowed among taxa in sampling rates. Essentially, this function allows taxa to be broken down into different possibly overlapping classes which have 'average' parameter values that are then combined to calcualte per-taxon parameters. For example, perhaps I think that taxa that live in a particular environment have a different characteristic sampling rate/probability, taxa of several different major clades have different characteristic sampling parameters and that there may be several temporal shifts in the characteristic extinction rate or sampling parameters. The classification IDs for the first two can be included as grp1 and grp2 and the hypothesized number of temporal breaks can be included as the n_tbins argument. A model where taxa differ in parameters across time, clades and environments will then be fit and the AIC calculated, so that it can be compared to other models. 

By default, the simple model where all taxa belong to a single class, with a single characteristic extinction rate and a single characteristic sampling parameter, is fit to the range data.

The timebins option will always allow for timebins with free-floating boundaries that are not defined a priori. The boundaries between time bins with different characteristic parameters will thus be additional parameters included in the AIC calculation. If you have the prior inclination that sampling/extinction changed at a particular point in time, then seperate the taxa that originated before and after that point as two different groups and include those classifications as a grp in the arguments.

This function does not implement the finite window of observation modification for dealing with data that leads up to the recent. This is planned for a future update, however. Thus, data input into this should be for taxa that have already gone extinct by the present and are not presently extant.

The timeData should be non-overlapping sequential intervals of roughly equal length. They should be in relative time, so the earliest interval should be 1 and it should increase as the intervals go up with age. This is so differences in interval numbers represents the same rough difference in interval timing. For example, taxa with listed as from North American Mammal assemblage zones, microfossil faunal zones or graptolite biozones can be given as long as they are correctly numbered in sequential. As a counter example, taxa resolved only to the Mesozoic and taxa resolved to the Turinean should not be included in the same input. I reccommend droping the taxa from less poorly resolved intervals.

The optimization of the maximum-likelihood solution for sampling probability is computationally slow because completeness (Pp in Foote, 1997) must be calculated for every parameter combination to obtain the likelihood, and that operation requires summation over a discrete probability distribution. To deal with this, by default this program will ignore the tails of the distributions, which decreases run time but also decreases accuracy of the likelihood calculation. For very flat likelihood surfaces, the accuracy of the likelihood calculation can be increased by increasing the threshold (I recommend a maximum threshold of 4). The default threshold is 1.

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

If the multi-class models are using, the element $pars will not be present, but there will be several different elements that sum the characteristic parameter components for each class. As noted in the $title, these should not be interpretated as the actual rates/probabilities of any real taxa but rather as components which must be assessed in combination with other classes to be meaningful. For example, for taxa of a given group in a given time bin, their extinction rate is the extinction rate component of that time bin times the extinction rate component of their group.

Completeness estimates are only output when model classes are not overlapping (and thus 'meaningful').
}
\references{
Foote, M. 1997. Estimating Taxonomic Durations and Preservation Probability. Paleobiology 23(3):278-300.

Foote, M., and D. M. Raup. 1996. Fossil preservation and the stratigraphic ranges of taxa. Paleobiology 22(2):121-140.
}
\author{David W. Bapst}
\seealso{
\code{\link{getSampRateCont}},\code{\link{sProb2sRate}},
}
\examples{
#Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
#we can do something very similar for the discrete time data (can be a bit slow)
SPres1<-getSampProbDisc(rangesDisc)
sProb<-SPres1[[2]][2]
print(sProb)
#note that that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we want the rate for srcTimePaleoPhy()
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
sRate<-sProb2sRate(sProb,int.length=1)
print(sRate)
#estimates that r=0.3... kind of low... (simulated sampling rate is 0.5)
\dontrun{
#Note: for real data, you may need to use an average int.length (no constant length)
#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins, AIC should be higher (less informative models)
randomgroup<-sample(1:2,nrow(rangesDisc[[2]]),replace=TRUE)
SPres2<-getSampProbDisc(rangesDisc,grp1=randomgroup)
SPres3<-getSampProbDisc(rangesDisc,n_tbins=2)
print(c(SPres1$AICc,SPres2$AICc,SPres3$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)
}
}
\keyword{datagen}
