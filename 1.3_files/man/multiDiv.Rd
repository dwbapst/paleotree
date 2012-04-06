\name{multiDiv}
\alias{multiDiv}
\alias{plotMultiDiv}
\title{Calculating Diversity Curves Across Multiple Datasets}
\description{Given a list of data, such as taxic range and phylogenetic trees, diversity curves are calculated for the same intervals across all the individual datasets. A median curve with 95 percent quantile bounds can be calculated and plotted.
}
\usage{
multiDiv(data, int.length = 1, plot = TRUE, split.int = TRUE, drop.ZLB=TRUE, plotLogRich = FALSE, timelims = NULL)
plotMultiDiv(results, plotLogRich = FALSE, timelims = NULL)
}
\arguments{
  \item{data}{A list where each element is a dataset formatted to be input in one of the diversity curve functions in this package}
  \item{int.length}{Time interval length}
  \item{plot}{Should the median diversity curve be plotted}
  \item{drop.ZLB}{Should terminal zero-length branches be dropped on phylogenetic datasets?}
  \item{split.int}{For discrete time data, should calculated/input intervals be split at discrete time interval boundaries? If FALSE, can create apparent artifacts in calculating the diversity curve.}
  \item{results}{The output of a previous run of multiDiv}
  \item{plotLogRich}{If true, taxic diversity plotted on log scale}
  \item{timelims}{Limits for the x (time) axis for diversity curve plots. Only affects plotting. Given as either NULL (the default) or as a vector of length two as for 'xlim' in the basic R function plot.}
}
\details{
This function is essentially a wrapper for the individual diversity curve functions included in this package (paleotree). It will intuitively decide whether datasets it is given are continuous-time taxic ranges, discrete-time (binned interval) taxic ranges or phylogenetic trees, as long as they are formatted as the individual functions request. A list that contains a mix of data types is entirely acceptable. A list of matrices output from simFossilTaxa is allowable, and considered as input for taxicDivCont. Data of an unknown type gives back an error.

The argument split.int splits intervals, if and only if discrete interval time data is included among the datasets. See the help file for taxicDivDisc to see an explanation of why split.int=TRUE by default is probably a good thing.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

The 'averaged' curve is actually the median rather than the mean as diversity counts are often highly skewed (in my experience).

The shaded certainty region around the median curve is the two-tailed 95 percent lower and upper quantiles, calculated from the observed data. It is not a true confidence interval, as it has no relaitonship to the standard error.
}
\value{
A list composed of three elements will be invisibly returned:
\item{int.times}{A two column matrix giving interval start and end times}
\item{div}{A matrix of measured diversities in particular intervals by rows, with each column representing a different dataset included in the input}
\item{median.curve}{A three column matrix, where the first column is the calculated median curve and the second and third columns are the 95 percent quantile upper and lower bounds}
}
\author{David W. Bapst}
\seealso{
The diversity curve functions used: \code{\link{phyloDiv}},\code{\link{taxicDivCont}},\code{\link{taxicDivDisc}},

Also see the function LTT.average.root() in the package TreeSim, which calculates an average LTT curve for multiple phylogenies.
}
\examples{
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
rangesCont<-sampleRanges(taxa,r=0.5)
rangesDisc<-binTimeData(rangesCont,int.length=1)
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#using multiDiv with very different data types
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=FALSE)
input<-list(rangesCont,rangesDisc,ttree)
multiDiv(input,plot=TRUE)

#using multiDiv with samples of trees
ttrees<-timePaleoPhy(cladogram,rangesCont,type="basic",randres=TRUE,ntrees=10,add.term=TRUE)
multiDiv(ttrees)
#uncertainty in diversity history is solely due to 
 #the random resolution of polytomies

#multiDiv can also take output from simFossilTaxa
#what do many simulations run under some conditions 'look' like on average?
set.seed(444)
taxa<-simFossilTaxa(p=0.3,q=0.1,nruns=20,maxtime=20,maxtaxa=100,plot=TRUE,min.cond=FALSE)
multiDiv(taxa)
#increasing cone of diversity! Even better on a log scale:
multiDiv(taxa,plotLogRich=TRUE)
}