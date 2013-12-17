\name{multiDiv}
\alias{multiDiv}
\alias{plotMultiDiv}
\title{Calculating Diversity Curves Across Multiple Datasets}
\description{Calculates multiple diversity curves from a list of datasets of taxic ranges and/or phylogenetic trees, for the same intervals, for all the individual datasets. A median curve with 95 percent quantile bounds is also calculated and plotted for each interval.
}
\usage{
multiDiv(data, int.length = 1, plot = TRUE, split.int = TRUE, drop.ZLB=TRUE, 
    drop.cryptic = FALSE, extant.adjust = 0.01, plotLogRich = FALSE, 
    timelims = NULL, int.times = NULL, plotMultCurves = FALSE,
    multRainbow = TRUE, divPalette = NULL)

plotMultiDiv(results, plotLogRich = FALSE, timelims = NULL, 
    plotMultCurves = FALSE, multRainbow = TRUE, divPalette = NULL)
}
\arguments{
  \item{data}{A list where each element is a dataset formatted to be input in one of the diversity curve functions in this package}
  \item{int.length}{Time interval length}
  \item{plot}{Should the median diversity curve be plotted}
  \item{drop.ZLB}{Should terminal zero-length branches be dropped on phylogenetic datasets?}
  \item{drop.cryptic}{If true, cryptic taxa are merged to form one taxon for estimating taxon curves. Only works for objects from simFossilTaxa.}
  \item{extant.adjust}{Amount of time to be added to extend start time for (0,0) bins for extant taxa}
  \item{split.int}{For discrete time data, should calculated/input intervals be split at discrete time interval boundaries? If FALSE, can create apparent artifacts in calculating the diversity curve.}
  \item{int.times}{An optional two-column matrix of the interval start and end times for calculating the diversity curve. If NULL, calculated internally. If given, the argument split.int and int.length are ignored.}
  \item{results}{The output of a previous run of multiDiv}
  \item{plotLogRich}{If true, taxic diversity plotted on log scale}
  \item{timelims}{Limits for the x (time) axis for diversity curve plots. Only affects plotting. Given as either NULL (the default) or as a vector of length two as for 'xlim' in the basic R function plot.}
  \item{plotMultCurves}{If TRUE, each individual diversity curve is plotted rather than the median diversity curve and 95 percent quantiles. FALSE by default.}
  \item{multRainbow}{If TRUE and plotMultCurves are both TRUE, each line is plotted as a different, randomized color using the function 'rainbow'. If FALSE, each line is plotted as a black line. This argument is ignored if divPalette is supplied.}
  \item{divPalette}{Can be used so users can pass a vector of chosen color identifiers for each diversity curve in 'data' which will take precedence over multRainbow. Must be the same length as the number of diversity curves supplied.}
}
\details{
This function is essentially a wrapper for the individual diversity curve functions included in paleotree. multiDiv will intuitively decide whether input datasets are continuous-time taxic ranges, discrete-time (binned interval) taxic ranges or phylogenetic trees, as long as they are formatted as required by the respective diversity curve functions. A list that contains a mix of data types is entirely acceptable. A list of matrices output by simFossilTaxa is allowable, and treated as input for taxicDivCont. Data of an unknown type gives back an error.

The argument split.int splits intervals, if and only if discrete interval time data is included among the datasets. See the help file for taxicDivDisc to see an explanation of why split.int=TRUE by default is probably a good thing.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

The 'averaged' curve is actually the median rather than the mean as diversity counts are often highly skewed (in this author's experience).

The shaded certainty region around the median curve is the two-tailed 95 percent lower and upper quantiles, calculated from the observed data. It is not a true probabilisitc confidence interval, as it has no relationship to the standard error.
}
\value{
A list composed of three elements will be invisibly returned:
\item{int.times}{A two column matrix giving interval start and end times}
\item{div}{A matrix of measured diversities in particular intervals by rows, with each column representing a different dataset included in the input}
\item{median.curve}{A three column matrix, where the first column is the calculated median curve and the second and third columns are the 95 percent quantile upper and lower bounds}
}
\author{David W. Bapst}
\seealso{
The diversity curve functions used include: \code{\link{phyloDiv}}, \code{\link{taxicDivCont}} and \code{\link{taxicDivDisc}}.

Also see the function LTT.average.root in the package TreeSim, which calculates an average LTT curve for multiple phylogenies, the functions mltt.plot in ape and ltt in phytools.
}
\examples{
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
rangesCont <- sampleRanges(taxa,r=0.5)
rangesDisc <- binTimeData(rangesCont,int.length=1)
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#using multiDiv with very different data types
ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=FALSE)
input <- list(rangesCont,rangesDisc,ttree)
multiDiv(input,plot=TRUE)

#using fixed interval times
multiDiv(input,int.times=rangesDisc[[1]],plot=TRUE)

#using multiDiv with samples of trees
ttrees <- timePaleoPhy(cladogram,rangesCont,type="basic",randres=TRUE,ntrees=10,add.term=TRUE)
multiDiv(ttrees)
#uncertainty in diversity history is solely due to 
   #the random resolution of polytomies

#multiDiv can also take output from simFossilTaxa
#what do many simulations run under some conditions 'look' like on average?
set.seed(444)
taxa <- simFossilTaxa(p=0.3,q=0.1,nruns=20,maxtime=20,maxtaxa=100,plot=TRUE,min.cond=FALSE)
multiDiv(taxa)
#increasing cone of diversity! 
#Even better on a log scale:
multiDiv(taxa,plotLogRich=TRUE)

#pure-birth example with simFossilTaxa
#note that conditioning is tricky
taxa <- simFossilTaxa(p=0.1,q=0,mintime=10,mintaxa=10,maxtime=50,maxtaxa=50,
    nruns=10,plot=TRUE)
multiDiv(taxa,plotLogRich=TRUE)

#compare many discrete diversity curves
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=20,maxtime=20,
    mintaxa=10,maxtaxa=100,plot=FALSE,min.cond=FALSE)
multiDiv(lapply(taxa,function(x) binTimeData(sampleRanges(x,r=0.5,min.taxa=1),int.length=1)))

layout(1)
}