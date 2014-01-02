\name{binTimeData}
\alias{binTimeData}
\title{Bin Simulated Temporal Ranges in Discrete Intervals}
\usage{
binTimeData(timeData, int.length = 1, start = NA, int.times = NULL)
}
\arguments{
  \item{timeData}{Two-column matrix of simulated first and
  last occurrences in absolute continuous time}

  \item{int.length}{Time interval length, default is 1 time
  unit}

  \item{start}{Starting time for calculating the
  intervals.}

  \item{int.times}{A two column matrix with the start and
  end times of the intervals to be used.}
}
\value{
A list containing: \item{int.times}{A 2 column matrix with
the start and end times of the intervals used; time
decreases relative to the present.} \item{taxon.times}{A 2
column matrix with the first and last occurrences of taxa
in the intervals listed in int.times, with numbers
referring to the row of int.times.}
}
\description{
Converts a matrix of simulated continuous-time first
occurrences and last occurrences for fossil taxa into first
and last occurrences given in some set of simulated/placed
discrete-time intervals, which is output along with
information of the dates of the given intervals.
}
\details{
This function takes a simulated matrix of per-taxon first
and last occurrences and, by dividing the time-scale into
time intervals of non-zero length, lists taxon occurrences
within those interval. By default, a set of sequential
non-overlapping time-interval of equal non-zero length are
used, with the length controlled by the argument
int.length.

Alternatively, a two column matrix of interval start and
end times to be used can be input via the argument
int.times. None of these intervals can have a duration
(temporal length) greater than zero. If a first or last
appearance in the input range data could fit into multiple
intervals (i.e. the input discrete time intervals are
overlapping), then the appearance data is placed in the
interval of the shortest duration. When output, the
interval times matrix (see below) will be sorted from first
to last.

This function is SPECIFICALLY for simulating the effect of
having a discrete timescale for analyses using simulations.
This function should not be used for non-simulations uses,
such as binning temporal occurrences for analyses of real
data. In those case, the temporal ranges (which, in real
data, will probably be given as discrete time intervals)
should already be tabulated within discrete intervals prior
to use in paleotree. The user should place the temporal
information in a list object, as described for the output
of binTimeData (see below).

As with many functions in the paleotree library, absolute
time is always decreasing, i.e. the present day is zero.
However, the numbering of intervals giving in the output
increases with time, as these are numbered relative to each
other, from first to last.

As of version 1.7 of paleotree, taxa which are extant as
indicated in timeData as being in a time interval bounded
(0,0), unless time-bins are preset using argument int.times
(prior to version 1.5 they were erroneously listed as NA).
}
\examples{
#Simulate some fossil ranges with simFossilTaxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont <- sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc <- binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
equalDiscInt <- taxicDivDisc(rangesDisc)

#example with pre-set intervals input (including overlapping)
presetIntervals <- cbind(c(1000,990,970,940),c(980,970,950,930))
rangesDisc1 <- binTimeData(rangesCont,int.times=presetIntervals)
taxicDivDisc(rangesDisc1)
#now let's plot diversity with (different) equal length intervals used above
taxicDivDisc(rangesDisc1,int.times=equalDiscInt[,1:2])

#example with extant taxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont <- sampleRanges(taxa,r=0.5,,modern.samp.prob=1)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc <- binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)

#example with pre-set intervals input (including overlapping)
presetIntervals <- cbind(c(40,30,20,10),c(30,20,10,0))
rangesDisc1 <- binTimeData(rangesCont,int.times=presetIntervals)
taxicDivDisc(rangesDisc1)
}
\seealso{
\code{\link{simFossilTaxa}}, \code{\link{sampleRanges}},
\code{\link{taxicDivCont}}
}

