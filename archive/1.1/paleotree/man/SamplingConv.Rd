\name{SamplingConv}
\alias{sProb2sRate}
\alias{sRate2sProb}
\alias{pqsRate2sProb}
\alias{qsProb2Comp}
\alias{qsRate2Comp}
\alias{probAnc}
\title{Converting Sampling Estimates}
\description{Various functions for converting between estimates of sampling in the fossil record}
\usage{
sProb2sRate(R, int.length = 1)
sRate2sProb(r, int.length = 1)
pqsRate2sProb(r, p, q, int.length = 1)
qsProb2Comp(R, q)
qsRate2Comp(r, q)
probAnc(p, q, R)
}
\arguments{
  \item{R}{Per-interval probability of sampling a taxon at least once}
  \item{r}{Instantaneous rate of sampling}
  \item{p}{Instantaneous rate of speciation (lambda)}
  \item{q}{Instantaneous rate of extinction (mu)}
  \item{int.length}{Length of Time Intervals}
}

\details{
This is a family of functions which all convert from some estimate of sampling to another estimate of sampling. Some of these also require estimates of the diversification rate. Diversification rates should always be the instanteous rate (often called the per-capita rates by paleontologists).

As with many models used in the paleotree library, it is generally assumed that the fossil record of interest is composed of discrete relatively-static taxonomic units which diversify mainly by budding cladogenesis, and that sampling is a very rare event well approximated by a model of exponential waiting times. The veracity of those assumptions is difficult to test and the senstivity of these analyses to relaxing those assumptions probably varies.

sProb2sRate() and sRate2sProb() give very rough conversions for converting between the probability of sampling once per time interval (R or "sProb" in this package as used in the references below) and the instantaneous rate of sampling per lineage/time unit ("sRate" or r). If you have estimates of the speciation and extinction rate, use pqsRate2sProb() instead for a more accurate estimate of R.

qsProb2Comp() and qsRate2Comp() are different calculations for "Pp" or the probability/proportion of taxa sampled in a clade. Theoretically, one could use it to extrapolate out the 'true' diversity, assuming the sampling rate model was correct. See Foote and Raup, 1996.

probAnc() obtains the probability of sampling indirect ancestors for a group as given by Foote (1996), assuming only budding cladogenesis, given the sampling probability and estimates of the diversification rate.

See the references below for a more detailed explanation of the methods and formulae used. The relevant equations are generally found in the appendices of those papers.
}
\value{
The converted sampling estimate, depending on the function used. See details above.
}
\references{
Foote, M. 1996. On the Probability of Ancestors in the Fossil Record. Paleobiology 22(2):141-151.

Foote, M. 1997. Estimating Taxonomic Durations and Preservation Probability. Paleobiology 23(3):278-300.

Foote, M. 2000. Origination and extinction components of taxonomic diversity: general problems. Pp. 74-102. In D. H. Erwin, and S. L. Wing, eds. Deep Time: Paleobiology's Perspective. The Paleontological Society, Lawrence, Kansas.

Foote, M., and D. M. Raup. 1996. Fossil preservation and the stratigraphic ranges of taxa. Paleobiology 22(2):121-140.

Solow, A. R., and W. Smith. 1997. On Fossil Preservation and the Stratigraphic Ranges of Taxa. Paleobiology 23(3):271-277.
}
\author{David W. Bapst, with advice from Michael Foote}
\seealso{
\code{\link{sampleRanges}},\code{\link{getSampRateCont}},\code{\link{getSampProbDisc}}
}
\examples{
sRate2sProb(r=0.5)
sProb2sRate(R=0.1)
pqsRate2sProb(r=0.5,p=0.1,q=0.1)
qsProb2Comp(R=0.1,q=0.1)
qsRate2Comp(r=0.1,q=0.1)
probAnc(p=0.1,q=0.1,R=0.5)
}
\keyword{datagen}
