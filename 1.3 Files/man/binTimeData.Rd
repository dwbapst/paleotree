\name{binTimeData}
\alias{binTimeData}
\title{Bin Temporal Ranges in Discrete Intervals}
\description{
Converts a matrix of continuous-time first occurances and last occurances for fossil taxa into first and last occurances given in discrete-time intervals, along with information the dates of the given intervals.
}
\usage{
binTimeData(timeData, int.length = 1, start = NA)
}
\arguments{
\item{timeData}{Two-column matrix of first and last occurrances in absolute continous time}
\item{int.length}{Time Interval Length, default is 1 time unit}
\item{start}{The starting time for calculating the intervals.}
}
\details{
This function takes a matrix of per-taxon first and last occurrances and divides those occurances into sequential non-overlappling time intervals, with each interval as long as int.length.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

The numbering of intervals given in the output increase with time, however, as these are in relative time.
}
\value{
a list containing:
\item{int.times}{A 2 column matrix with the start and end times of the intervals used.}
\item{taxon.times}{A 2 column matrix with the first and last occurances of taxa in the intervals used, with numbers referring to the row of int.times.}
}
\author{David W. Bapst}
\seealso{
\code{\link{simFossilTaxa}},\code{\link{sampleRanges}},\code{\link{taxicDivCont}}
}
\examples{
dev.new()
#Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)
}
