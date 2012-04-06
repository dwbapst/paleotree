\name{getSampProbDisc}
\alias{getSampProbDisc}
\title{Fit Models of Sampling Probability to Discrete-Interval Taxon Ranges}
\description{Uses ML to find the best-fit parameters for models of sampling probability and extinction rate, given a set of discrete-interval taxon ranges from the fossil record. This function can fit models where there are different groupings of taxa with different parameters and different free-floating time intervals with different parameters.}
\usage{
getSampProbDisc(timeData, n_tbins = 1, grp1 = NA, grp2 = NA, est_only = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timeData}{A 2 column matrix with the first and last occurances of taxa given in relative time intervals. If a list of length two is given for timeData, such as would be expected if the output of binTimeData was directly input, the second element is used.}
  \item{n_tbins}{Number of time bins with different sampling/extinction parameters}
  \item{grp1}{A vector of integers or characters, the same length as the number of taxa in timeData, where each taxon-wise element gives the group ID of the taxon for the respective row of timeData}
  \item{grp2}{A vector of integers or characters, the same length as the number of taxa in timeData, where each taxon-wise element gives the group ID of the taxon for the respective row of timeData}
  \item{est_only}{If true, function will give back a matrix of ML extinction rates and sampling probabilities per species rather than usual output (see below)}
}
\details{
This function uses maximum-likelihood solutions found by Foote (1997). These analyses are ideally applied to data from single stratigraphic section but can potentially be applicable to regional or global datasets, although the behavior of those datasets is less well understood.

getSampProbDisc allows for a considerable level of versatility in terms of the degree of variation allowed among taxa in sampling rates. Essentially, this function allows taxa to be broken down into different possibly overlapping classes which have 'average' parameter values that are then combined to calcualte per-taxon parameters. For example, perhaps I think that taxa that live in a particular environment have a different characteristic sampling rate/probability, taxa of several different major clades have different characteristic sampling parameters and that there may be several temporal shifts in the characteristic extinction rate or sampling parameters. The classification IDs for the first two can be included  as per-taxon vectors (of either character or integers) as grp1 and grp2 and the hypothesized number of temporal breaks can be included as the n_tbins argument. A model where taxa differ in parameters across time, clades and environments will then be fit and the AIC calculated, so that it can be compared to other models. 

By default, the simple model where all taxa belong to a single class, with a single characteristic extinction rate and a single characteristic sampling parameter, is fit to the range data.

The timebins option allows for timebins with free-floating boundaries that are not defined a priori. The boundaries between time bins with different characteristic parameters will thus be additional parameters included in the AIC calculation. If you have the prior inclination that sampling/extinction changed at a particular point in time, then seperate the taxa that originated before and after that point as two different groups and include those classifications as a grp in the arguments.

This function does not implement the finite window of observation modification for dealing with data that leads up to the recent. This is planned for a future update, however. For now, data input into this function should be for taxa that have already gone extinct by the modern and are not presently extant.

The timeData should be non-overlapping sequential intervals of roughly equal length. They should be in relative time, so the earliest interval should be 1 and it should increase as the intervals go up with age. This is so differences in interval numbers represents the same rough difference in interval timing. For example, a dataset where all taxa are listed from a set of sequential intervals of similar length, such as North American Mammal assemblage zones, microfossil faunal zones or graptolite biozones can be given as long as they are correctly numbered in sequential order in the input. As a counter example, a dataset which includes taxa resolved only to intervals as wide as the whole Jurassic and taxa resolved to biozones within the Jurassic should not be included in the same input. Drop the taxa from less poorly resolved intervals from such datasets if you want to apply this function, as long as this retains a large enough sample of taxa from the sequential intervals. Note that taxicDivDisc and the "bin_" timescaling methods do not require that intervals be truly sequential (they can be overlapping; see their helpfiles). The output from binTimeData is always sequential, at least by default.

Please check the $message element of the output to make sure that convergence occurred.
}
\value{
If est_only=T, a matrix of per-taxon sampling and extinction parameters is output.

If est_only=F (default), then the output is a list:

\item{Title}{Gives details of the analysis, such as the number and type of parameters included and the number of taxa analyzed}
\item{pars}{Maximum-likelihood parameters of the sampling model, per class of taxa fit}
\item{SMax}{The maximum support (log-likelihood) value}
\item{AICc}{The second-order Akaike's Information Criterion, corrected for small sample sizes}
\item{message}{Messages output by optim(); check to make sure that model convergence occurred}

If the multi-class models are using, the element $pars will not be present, but there will be several different elements that sum the characteristic parameter components for each class. As noted in the $title, these should not be interpretated as the actual rates/probabilities of any real taxa but rather as components which must be assessed in combination with other classes to be meaningful. For example, for taxa of a given group in a given time bin, their extinction rate is the extinction rate component of that time bin times the extinction rate component of their group. Completeness estimates (Pp) will be output with these parameters as long as classes are not overlapping, as those estimates would not otherwise refer to meaningful groups of taxa.
}
\references{
Foote, M. 1997. Estimating Taxonomic Durations and Preservation Probability. Paleobiology 23(3):278-300.

Foote, M., and D. M. Raup. 1996. Fossil preservation and the stratigraphic ranges of taxa. Paleobiology 22(2):121-140.
}
\author{David W. Bapst, with considerable advice from Michael Foote.}
\seealso{
\code{\link{getSampRateCont}},\code{\link{sProb2sRate}},\code{\link{qsProb2Comp}},
}
\examples{
#Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#now, get an estimate of the sampling rate (we set it to 0.5 above)
#for discrete data we can estimate the sampling probability per interval (R)
	#i.e. this is not the same thing as the instantaneous sampling rate (r)
#can use sRate2sProb to see what we would expect
sRate2sProb(r=0.5)
#expect R = ~0.39
#now we can use maximum likelihood to taxon ranges to get sampling probability
SPres1<-getSampProbDisc(rangesDisc)
sProb<-SPres1[[2]][2]
print(sProb)
#est. R = ~0.42; not too off what we would expect!
#for the src based timescaling methods, we want an estimate of the instantaneous samp rate
#we can use sProb2sRate() to get the rate. We will also need to also tell it the int.length
sRate<-sProb2sRate(sProb,int.length=1)
print(sRate)
#estimates that r=0.54... Not bad!
#Note: for real data, you may need to use an average int.length (no constant length)

\dontrun{
#this data was simulated under homogenous sampling rates, extinction rates
#if we fit a model with random groups and allow for multiple timebins
	#AIC should be higher (less informative models)
randomgroup<-sample(1:2,nrow(rangesDisc[[2]]),replace=TRUE)
SPres2<-getSampProbDisc(rangesDisc,grp1=randomgroup)
SPres3<-getSampProbDisc(rangesDisc,n_tbins=2)
print(c(SPres1$AICc,SPres2$AICc,SPres3$AICc))
#and we can see the most simple model has the lowest AICc (most informative model)

#testing temporal change in sampling rate
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=100,maxtaxa=125,maxtime=1000,maxExtant=0,plot=T)
#let's see what the 'true' diversity curve looks like in this case
#simulate two sets of ranges at r=0.7 and r=0.1
rangesCont<-sampleRanges(taxa,r=1.1)
rangesCont2<-sampleRanges(taxa,r=0.2)
#now make it so that taxa which originated after 850 have r=0.1
rangesCont[taxa[,3]<850,]<-rangesCont2[taxa[,3]<850,]
rangesDisc<-binTimeData(rangesCont)
#lets plot the diversity curve
taxicDivDisc(rangesDisc)
SPres1<-getSampProbDisc(rangesDisc)
SPres2<-getSampProbDisc(rangesDisc,n_tbins=2)
#compare the AICc of the models
print(c(SPres1$AICc,SPres2$AICc)) #model 2 looks pretty good
#when does it find the break in time intervals?
print(rangesDisc[[1]][SPres2$t_ends[2],1])
#not so great: estimates 940, not 850 
	#but look at the diversity curve: most richness in bin 1 is before 940
	#might have found the right break time otherwise...
#the parameter values it found are less great. Finds variation in q	
}
}
\keyword{datagen}
