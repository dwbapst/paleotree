\name{getSampRateCont}
\alias{getSampRateCont}
\title{Fit Models of Sampling Rates to Continuous-Time Taxon Ranges}
\description{Uses maximum likelihood to estimate instantaneous sampling and extinction rates, given a set of continuous-time taxon ranges from the fossil record
}
\usage{
getSampRateCont(timeData, n_tbins = 1, grp1 = NA, grp2 = NA, threshold = 0.1, est_only = FALSE, iterations = 1000000, initial = 0.5)
}
\arguments{
  \item{timeData}{Two-column matrix of per-taxon first and last occurrances in absolute continous time}
  \item{n_tbins}{Number of time bins with different sampling/extinction parameters}
  \item{grp1}{A vector of integers or characters, the same length as the number of taxa in timeData, where each taxon-wise element gives the group ID of the taxon for the respective row of timeData}
  \item{grp2}{A vector of integers or characters, the same length as the number of taxa in timeData, where each taxon-wise element gives the group ID of the taxon for the respective row of timeData}
  \item{threshold}{The smallest allowable range (measured difference in the FAD and LAD of a taxon). Ranges below this size will be treated as "one-hit" sampling events.}
  \item{est_only}{If true, function will give back a matrix of ML extinction rates and sampling probabilities per species rather than usual output (see below)}
  \item{iterations}{Maximum number of iterations the optimizer is run for}
  \item{initial}{Values used as initial parameter value to initiate the optimizing function. The same starting value is used for all parameters}
}
\details{
This function uses maximum-likelihood solutions obtained by Foote (1997). These analyses are ideally applied to data from single stratigraphic section but potentially are applicable to regional or global datasets, although the behavior of those datasets is less well understood.

getSampRateCont allows for a considerable level of versatility in terms of the variation allowed among taxa in model parameters (extinction rate and sampling probability/rate). Essentially, taxa are divided into different (possibly overlapping) classes which have 'average' parameter values. These average parameters are multiplicatively-combined to calcualte per-taxon parameters. For example, perhaps a user hypotheses taxa that live in a particular environment have different characteristic sampling rates/probabilities, that taxa of several different major clades have different characteristic sampling parameters and that there may be several temporal shifts in the characteristic extinction rate. The classification IDs for the first two can be included  as per-taxon vectors (of either character or integers) as grp1 and grp2 and the hypothesized number of temporal breaks included as the n_tbins argument. A model where taxa differ in parameters across time, clades and environments will thus be fit and the AIC calculated, such that this model's fit can be compared to other (probably simpler) models. 

By default, a simple model where all taxa belong to a single class, with a single characteristic extinction rate and a single characteristic sampling parameter, is fit to the range data.

The timebins option allows for time intervals with different characteristic model parameters which are not defined a priori. The boundaries between time intervals with different characteristic parameters will thus be additional free-floating parameters included in the AIC calculation. If you have the prior inclination that sampling/extinction changed at a specific point in time, then seperate the taxa that originated before and after that point as two different groups and include those classifications as a group in the arguments.

This function does not implement the finite window of observation modification for dealing with data that leads up to the recent (Foote, 1997). This is planned for a future update, however. For now, data input into this function should be for taxa that have already gone extinct by the modern and are not presently extant.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

Please check the $message and $convergence elements of the output to make sure that convergence occurred. The likelihood surface can be very flat in some cases, particularly for small datasets (<100 taxa). If convergence was not reached, a warning message is communicated. If the optimizer does not converge, consider increasing iterations or changing the starting values.
}
\value{
If est_only=T, a matrix of per-taxon sampling and extinction parameters is output.

If est_only=F (default), then the output is a list:

\item{Title}{Gives details of the analysis, such as the number and type of parameters included and the number of taxa analyzed}
\item{parameters}{Maximum-likelihood parameters of the sampling model}
\item{log.likelihood}{The maximum support (log-likelihood) value}
\item{AICc}{The second-order Akaike's Information Criterion, corrected for small sample sizes}
\item{convergence}{A number indicating status of convergence; if 0 then convergence was reached; see help file for optim for the respective error list}
\item{message}{Messages output by optim(); check to make sure that model convergence occurred}

If multi-class models are fit, the element $parameters will not be present, instead there will be several different elements that simple describe the characteristic parameter 'components' for each class, rather than representing the parameters of actual taxa in that class. As noted in the $title, these should not be interpretated as the actual rates/probabilities of any real taxa but rather as factors which must be multiplied in combination with the estimates for other classes to be meaningful. For example, for taxa of a given group in a given time bin, their extinction rate is the extinction rate component of that time bin times the extinction rate component of their group. Completeness estimates will be output with these parameters as long as classes are not overlapping, as those estimates would not otherwise refer to meaningful groups of taxa.
}
\references{
Foote, M. 1997 Estimating Taxonomic Durations and Preservation Probability. \emph{Paleobiology} \bold{23}(3):278--300.

Foote, M., and D. M. Raup. 1996 Fossil preservation and the stratigraphic ranges of taxa. \emph{Paleobiology} \bold{22}(2):121--140.
}
\author{David W. Bapst}
\seealso{
\code{\link{getSampProbDisc}}, \code{\link{sRate2sProb}}, \code{\link{qsRate2Comp}}
}
\examples{
#Simulate some fossil ranges with simFossilTaxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges
rangesCont <- sampleRanges(taxa,r=0.5)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
SRres1 <- getSampRateCont(rangesCont)
print(SRres1)   #that's all the results...

#pulling out the sampling rate from the results
sRate <- SRres1[[2]][2]
print(sRate)	#estimates that sRate=~0.4 (not too bad...)

#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins
	#AIC should be higher (less informative)
randomgroup <- sample(1:2,nrow(rangesCont),replace=TRUE)
SRres2 <- getSampRateCont(rangesCont,grp1=randomgroup)
SRres3 <- getSampRateCont(rangesCont,n_tbins=2)
SRres4 <- getSampRateCont(rangesCont,n_tbins=3,grp1=randomgroup)
print(c(SRres1$AICc,SRres2$AICc,SRres3$AICc,SRres4$AICc))
#the most simple model (the first value) has the lowest AICc (most informative!)
}