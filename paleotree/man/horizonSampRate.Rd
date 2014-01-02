\name{horizonSampRate}
\alias{horizonSampRate}
\title{Estimate Sampling Rate from Sampling Horizon Data (Solow and Smith, 1997)}
\usage{
horizonSampRate(sampOcc = NULL, durations = NULL, nCollections = NULL)
}
\arguments{
  \item{SampOcc}{A list with the number of elements equal
  to the number of taxa, and each element of the list being
  a numerical vector with the length equal to the number of
  collections for each taxon, and each value equal to the
  precise date of that fossil's time of collection. These
  dates do not need to be ordered. If not supplied, the
  elements durations and nCollections must be supplied.}

  \item{durations}{A vector of precise durations in
  continuous time, with the length equal to the number of
  taxa. If not suppliedthis is calculated from SampOcc,
  which must be supplied.}

  \item{nCollections}{A vector of integers representing the
  number of collections for each taxon in the input
  durations. If not supplied this is calculated from
  SampOcc, which must be supplied.}
}
\value{
Returns the instantaneous sampling (in lineage*time-units)
as a single numerical value. Note that this is the
instantaneous sampling rate and not the probability of
sampling a taxon per interval.
}
\description{
This function implements the exact maximum likelihood
estimator for the instantaneous sampling rate from Solow
and Smith (1997, Paleobiology), which is based on the
relationship between the number of collections for a set of
tax and their durations (known precisely in continuous
time).
}
\details{
Given a dataset of taxa with a vector N, representing the
number of collections for each taxon, and a vector D,
giving the precise duration for each taxon, we can use the
following maximum likelihood estimator from Solow and Smith
(1997) to obtain the instantaneous sampling rate:

\eqn{samplingRate = (sum(N-1)^2)/(sum(D)*sum(N))}

This method is exclusively for datasets with very precisely
dated horizons, such as microfossils from deep sea cores
with very precise age models. The first and last appearance
must be known very precisely to provide an equally precise
estimate of the duration. Most datasets are not precise
enough for this method, due to chronostratigraphic
uncertainty. However, note that the age of individual
collections other than the first and last appearance dates
do not need to be known: its only the number of collections
that matters.
}
\examples{
#can simulate this type of data with sampleRanges
    # just set ranges.only = FALSE
#let's try a simulation example:
set.seed(444)
taxa <- simFossilTaxa(p = 0.1,q = 0.1,nruns = 1,mintaxa = 20,maxtaxa = 30,maxtime = 1000,maxExtant = 0)
sampledOccurrences <- sampleRanges(taxa,r = 0.5,ranges.only = FALSE)

# now try with horizonSampRate
horizonSampRate(sampOcc = sampledOccurrences)

# but we could also try with the *other* inputs
   # useful because some datasets we may only have durations
   # and number of sampling events for
filtered <- sampledOccurrences[!is.na(sampledOccurrences)]
dur <- sapply(filtered,max) - sapply(filtered,min)
nCol <- sapply(filtered,length)
# supply as durations and nCollections
horizonSampRate(durations = dur, nCollections=nCol)
}
\references{
Solow, A. R., and W. Smith. 1997. On Fossil Preservation
and the Stratigraphic Ranges of Taxa. \emph{Paleobiology}
23(3):271-277.
}
\seealso{
Duration frequency methods (Foote and Raup, 1996; Foote,
1997) use ranges alone to estimate sampling parameters,
implemented in \code{\link{durationFreq}}.

Also see the conversion functions for sampling parameters
at \code{\link{SamplingConv}}.
}

