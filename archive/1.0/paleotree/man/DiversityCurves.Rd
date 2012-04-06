\name{DiversityCurves}
\alias{taxicDivCont}
\alias{taxicDivDisc}
\alias{phyloDiv}
\title{Diversity Curves}
\description{Plots diversity curves based on taxic range data (in both discrete and continuous time) and for phylogenies.}
\usage{
	taxicDivCont(timeData, int.length = 1, int.times = NULL, plot = T, output = F)
	taxicDivDisc(timeList, int.times = NULL, plot = T, output = F)
	phyloDiv(ttree, int.length = 1, int.times = NULL, output = F, plot = T, plotLogRich = F, drop.ZLB=T)
}
\arguments{
  \item{timeData}{Two-column matrix giving the per-taxon first and last appearances in absolute time}
  \item{timeList}{A list giving interval data. See details below.}
  \item{ttree}{a time-scaled phylogeny object of class phylo}
  \item{int.length}{The length of intervals used to make the diversity curve. Ignored if int.times is given.}
  \item{int.times}{A two-column matrix of the interval start and end times for calculating the diversity curve. If NULL, calculated internally.}
  \item{output}{If true, per-interval diversity data matrix is output, see details below}
  \item{plot}{If true, diversity curve is plotted}
  \item{plotLogRich}{If true, taxic diversity plotted on log scale}
  \item{drop.ZLB}{If true, zero-length terminal branches are dropped from the input tree}
}
\details{
Background: Diversity curves are just plots of how species richness changes over time for a particular group of organisms. In the fossil record, these are based on per-taxon range data while more recently molecular phylogenies have been used to plot lineage-through-time plots (LTTs). Neither of these are without their weaknesses; reconstructing the true history of biodiversity is exceedingly difficult for any group.

The diversity curves produced by these functions will always measure diversity within binned time intervals. For continuous-time data or phylogenies, one could decrease the int.length used to get what is essentially almost an 'instantaneous' estimate of diversity. This is warned against, however, as most historical diversity data will have some time-averaging or uncertain temporal resolution and thus should not be used to calculate instantaneous estimates of diversity.


As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

As diversity is counted within binned intervals, standing diversity may be somewhat lower than the measured quantity, particularly if intervals longer than the mean duration of taxa is used. This will be an issue with all diversity curve functions, but particularly the discrete-time variant.

The timeList object should be a list composed of two matrices, one giving by-interval start and end times (in absolute time), the other giving by-taxon first and last appearances in those intervals, numbered as the rows. Absolute time should be decreasing, while the intervals should be numbered so that the number increases with time. Unlike getSampProbDisc, the intervals can be overlapping. See the function of binTimeData for more information.

phyloDiv will resolve polytomies to be dichotomous nodes seperated by zero-length branches prior to calculating the diversity curve. There is not option to change this, but it should not affect the use of the function. phyloDiv will also drop zero-length terminal branches, as with the function dropZLB. This is done by default but can be turned off by setting the argument drop.zlb to false.
}
\value{
If output=TRUE, then these functions will output a three-column matrix, where the first two columns are interval start and end times and the third column is the number of taxa/lineages counted in that interval.
}
\author{David W. Bapst}
\seealso{
\code{\link{multiDiv}},\code{\link{ltt.plot}},

Also see the function ltt() in the package phytools. The packages TreeSim and LASER also have similar functions.
}
\examples{
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa[,3:4])
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#plot the diversity curve based on the sampled ranges
layout(matrix(1:2,,2))
taxicDivCont(rangesCont)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#plot with taxicDivDisc()
taxicDivDisc(rangesDisc)
#compare to the continuous time diversity curve
#Now let's make a tree using taxa2phylo()
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
phyloDiv(tree)

#a simple example with phyloDiv
set.seed(444)
tree<-rtree(100)
phyloDiv(tree)

#a neat example of using phyDiv with timeSliceTree to simulate doing molecular-phylogeny studies of diverification... in the past
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
taxicDivCont(taxa[,3:4])
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=TRUE)
}
\keyword{datagen}