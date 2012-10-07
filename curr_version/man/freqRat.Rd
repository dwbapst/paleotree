\name{freqRat}
\alias{freqRat}
\title{Frequency Ratio Method for Estimating Sampling Probability}
\description{Estimate per-interval sampling probability in the fossil record from a set of discrete-interval taxon ranges using the frequency-ratio method described by Foote and Raup (1996).}
\usage{
freqRat(timeData, plot = FALSE)
}
\arguments{
  \item{timeData}{A 2 column matrix with the first and last occurances of taxa given in relative time intervals. If a list of length two is given for timeData, such as would be expected if the output of binTimeData was directly input, the second element is used.}
  \item{plot}{If true, the histogram of observed taxon ranges is plotted, with frequencies on a linear scale}
}
\details{
This function uses the frequency ratio ("freqRat") method of Foote and Raup (1996) to estimate the per-interval sampling rate for a set of taxa. This method assumes that intervals are of fairly similar length and that taxonomic extinction and sampling works similar to homogenous Poisson processes. These analyses are ideally applied to data from single stratigraphic section but potentially are applicable to regional or global datasets, although the behavior of those datasets is less well understood.

The frequency ratio is a simple relationship between the number of taxa observed only in a single time interval (also known as singletons), the number of taxa observed only in two time intervals and the number of taxa observed in three time intervals. These respective frequencies, respectively f1, f2 and f3 can then be related to the per-interval sampling probability with the following expression.

\eqn{Sampling.Probability = (f2^2)/(f1*f3)}

This frequency ratio is generally referred to as the 'freqRat' in paleobiological literature.

It is relatively easy to visually test if range data fits expectation that true taxon durations are exponentially distributed by plotting the frequencies of the observed ranges on a log scale: data beyond the 'singletons' category should have a linear slope, implying that durations were originally exponentially distributed. (Note, a linear scale is used for the plotting diagram made by this function when 'plot' is TRUE, so that plot cannot be used for this purpose.)

The accuracy of this method tends to be poor at small interval length and even relatively large sample sizes. A portion at the bottom of the examples in the help file examine this issue in greater detail with simulations. This package author recommends using the ML method developed in Foote (1997) instead, which is usable via the function \code{\link{getSampProbDisc}}.
}
\value{
This function returns the per-interval sampling probability as the "freqRat".
}
\references{
Foote, M. 1997 Estimating Taxonomic Durations and Preservation Probability. \emph{Paleobiology} \bold{23}(3):278--300.

Foote, M., and D. M. Raup. 1996 Fossil preservation and the stratigraphic ranges of taxa. \emph{Paleobiology} \bold{22}(2):121--140.
}
\author{David W. Bapst}
\seealso{
\code{\link{getSampProbDisc}}, \code{\link{getSampRateCont}}, \code{\link{sProb2sRate}}, \code{\link{qsProb2Comp}}
}
\examples{
#Simulate some fossil ranges with simFossilTaxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=100,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges
rangesCont <- sampleRanges(taxa,r=0.1)
#Now let's use binTimeData to bin in intervals of 5 time units
rangesDisc <- binTimeData(rangesCont,int.length=5)

#now, get an estimate of the sampling rate (we set it to 0.5 above)
#for discrete data we can estimate the sampling probability per interval (R)
    #i.e. this is not the same thing as the instantaneous sampling rate (r)
#can use sRate2sProb to see what we would expect
sRate2sProb(r=0.1,int.length=5)
#expect R = ~0.39

#now we can use maximum likelihood to taxon ranges to get sampling probability
SampProb <- freqRat(rangesDisc,plot=TRUE)
SampProb
#est. R = ~0.25


\dontrun{
#################
#The following example code (which is not run by default) examines how 
	#the freqRat estimates vary with sample size, interval length
	#and compare it to using getSampProbDisc

#how good is the freqRat at 20 sampled taxa on avg?
set.seed(444)
r<-runif(100)
int.length=1
R<-sapply(r,sRate2sProb,int.length=1)	#estimate R from r, assuming stuff like p=q
ntaxa<-freqRats<-numeric()
for(i in 1:length(r)){
	#assuming budding model
	taxa<-simFossilTaxa_SRCond(r=r[i],avgtaxa=20,p=0.1,q=0.1,maxExtant=0)
	ranges<-sampleRanges(taxa,r=r[i])
	timeList<-binTimeData(ranges,int.length=int.length)
	ntaxa[i]<-nrow(timeList[[2]])
	freqRats[i]<-get_freqRat(timeList)
	}
plot(R,freqRats);abline(0,1)
#without the gigantic artifacts bigger than 1...
plot(R,freqRats,ylim=c(0,1));abline(0,1)
#very worrisome lookin'!

#how good is it at 100 sampled taxa on average?
set.seed(444)
r<-runif(100)
int.length=1
R<-sapply(r,sRate2sProb,int.length=1)
ntaxa<-freqRats<-numeric()
for(i in 1:length(r)){
	#assuming budding model
	taxa<-simFossilTaxa_SRCond(r=r[i],avgtaxa=100,p=0.1,q=0.1,maxExtant=0)
	ranges<-sampleRanges(taxa,r=r[i])
	timeList<-binTimeData(ranges,int.length=int.length)
	ntaxa[i]<-nrow(timeList[[2]])
	freqRats[i]<-get_freqRat(timeList)
	}
plot(R,freqRats,ylim=c(0,1));abline(0,1)
#not so hot, eh?

################
#LETS CHANGE THE TIME BIN LENGTH!

#how good is it at 100 sampled taxa on average, with longer time bins?
set.seed(444)
r<-runif(100)
int.length<-10
R<-sapply(r,sRate2sProb,int.length=int.length)
ntaxa<-freqRats<-numeric()
for(i in 1:length(r)){
	#assuming budding model
	taxa<-simFossilTaxa_SRCond(r=r[i],avgtaxa=100,p=0.1,q=0.1,maxExtant=0)
	ranges<-sampleRanges(taxa,r=r[i])
	timeList<-binTimeData(ranges,int.length=int.length)
	ntaxa[i]<-nrow(timeList[[2]])
	freqRats[i]<-get_freqRat(timeList)
	}
plot(R,freqRats,ylim=c(0,1));abline(0,1)
#things get more accurate as interval length increases... odd, eh?

#how good is it at 20 sampled taxa on average, with longer time bins?
set.seed(444)
r<-runif(100)
int.length<-10
R<-sapply(r,sRate2sProb,int.length=int.length)
ntaxa<-freqRats<-numeric()
for(i in 1:length(r)){
	#assuming budding model
	taxa<-simFossilTaxa_SRCond(r=r[i],avgtaxa=20,p=0.1,q=0.1,maxExtant=0)
	ranges<-sampleRanges(taxa,r=r[i])
	timeList<-binTimeData(ranges,int.length=int.length)
	ntaxa[i]<-nrow(timeList[[2]])
	freqRats[i]<-get_freqRat(timeList)
	}
plot(R,freqRats,ylim=c(0,1));abline(0,1)
#still not so hot at low sample sizes, even with longer bins

########################
#ML METHOD

#how good is the ML method at 20 taxa, 1 time-unit bins?
set.seed(444)
r<-runif(100)
int.length<-1
R<-sapply(r,sRate2sProb,int.length=int.length)
ntaxa<-ML_sampProb<-numeric()
for(i in 1:length(r)){
	#assuming budding model
	taxa<-simFossilTaxa_SRCond(r=r[i],avgtaxa=20,p=0.1,q=0.1,maxExtant=0)
	ranges<-sampleRanges(taxa,r=r[i])
	timeList<-binTimeData(ranges,int.length=int.length)
	ntaxa[i]<-nrow(timeList[[2]])
	ML_sampProb[i]<-getSampProbDisc(timeList)[[2]][2]
	}
plot(R,ML_sampProb);abline(0,1)
#Not so great due to likelihood surface ridges, but it returns values between 0-1

#how good is the ML method at 100 taxa, 1 time-unit bins?
set.seed(444)
r<-runif(100)
int.length<-1
R<-sapply(r,sRate2sProb,int.length=int.length)
ntaxa<-ML_sampProb<-numeric()
for(i in 1:length(r)){
	#assuming budding model
	taxa<-simFossilTaxa_SRCond(r=r[i],avgtaxa=100,p=0.1,q=0.1,maxExtant=0)
	ranges<-sampleRanges(taxa,r=r[i])
	timeList<-binTimeData(ranges,int.length=int.length)
	ntaxa[i]<-nrow(timeList[[2]])
	ML_sampProb[i]<-getSampProbDisc(timeList)[[2]][2]
	}
plot(R,ML_sampProb);abline(0,1)
#Oh, fairly nice, although still a biased uptick as R gets larger

}
}
