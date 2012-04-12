\name{DiversityCurves}
\alias{taxicDivCont}
\alias{taxicDivDisc}
\alias{phyloDiv}
\title{Diversity Curves}
\description{Plots diversity curves based on taxic range data (in both discrete and continuous time) and for phylogenies.}
\usage{
taxicDivCont(timeData, int.length = 1, int.times = NULL, plot = TRUE, plotLogRich = FALSE, timelims = NULL, drop.cryptic = FALSE)

taxicDivDisc(timeList, int.times = NULL, plot = TRUE, split.ints=TRUE, plotLogRich = FALSE, timelims = NULL)

phyloDiv(tree, int.length = 1, int.times = NULL, plot = TRUE, plotLogRich = FALSE, drop.ZLB = TRUE, timelims = NULL)
}
\arguments{
  \item{timeData}{Two-column matrix giving the per-taxon first and last appearances in absolute time. The matrix resulting from a single simFossilTaxon run can also be supplied to taxicDivCont, in which case all columns except the third and the fourth are ignored.}
  \item{timeList}{A list giving interval data. See details below.}
  \item{tree}{a time-scaled phylogeny object of class phylo}
  \item{int.length}{The length of intervals used to make the diversity curve. Ignored if int.times is given.}
  \item{int.times}{A two-column matrix of the interval start and end times for calculating the diversity curve. If NULL, calculated internally.}
  \item{split.int}{For discrete time data, should calculated/input intervals be split at discrete time interval boundaries? If FALSE, can create apparent artifacts in calculating the diversity curve. See below.}
  \item{plot}{If true, diversity curve is plotted}
  \item{plotLogRich}{If true, taxic diversity plotted on log scale}
  \item{drop.cryptic}{If true, cryptic taxa are merged to form one taxon for estimating taxon curves. Only works for objects from simFossilTaxa.}
  \item{timelims}{Limits for the x (time) axis for diversity curve plots. Only affects plotting. Given as either NULL (the default) or as a vector of length two as for 'xlim' in the basic R function plot.}
  \item{drop.ZLB}{If true, zero-length terminal branches are dropped from the input tree}
}
\details{
Background: Diversity curves are plots of species/taxon/lineage richness over time for a particular group of organisms. For paleontological studies, these are generally based on per-taxon range data while more recently in evolutionary biology, molecular phylogenies have been used to calculate lineage-through-time plots (LTTs). Neither of these are without their weaknesses; reconstructing the true history of biodiversity is exceedingly difficult for any group.

The diversity curves produced by these functions will always measure diversity within binned time intervals (and plot them as rectangular bins). For continuous-time data or phylogenies, one could decrease the int.length used to get what is essentially an 'instantaneous' estimate of diversity. This is warned against, however, as most historical diversity data will have some time-averaging or uncertain temporal resolution and thus probably is not finely-resolved enough to calculate instantaneous estimates of diversity.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

As diversity is counted within binned intervals, the true standing diversity may be somewhat lower than the measured (observed) quantity, particularly if intervals are longer than the mean duration of taxa is used. This will be an issue with all diversity curve functions, but particularly the discrete-time variant. For diversity data in particularly large discrete time intervals, plotting this data in smaller bins which do not line up completely with the original intervals will create a 'spiky' diversity curve, as these smaller intersecting bins will have a large number of taxa which may have been present in either of the neighboring intervals. This will give these small bins an apparently high estimated standing diversity. This artifact is avoided with the default setting split.int=TRUE, which will split any input or calculated intervals so that they start and end at the boundaries of the discrete-time range bins.

The timeList object should be a list composed of two matrices, one giving by-interval start and end times (in absolute time), the other giving by-taxon first and last appearances in those intervals, numbered as the rows. Absolute time should be decreasing, while the intervals should be numbered so that the number increases with time. Unlike getSampProbDisc, the intervals can be overlapping. See the function of binTimeData for more information.

phyloDiv will resolve polytomies to be dichotomous nodes seperated by zero-length branches prior to calculating the diversity curve. There is not option to change this, but it should not affect the use of the function. phyloDiv will also drop zero-length terminal branches, as with the function dropZLB. This is done by default but can be turned off by setting the argument drop.zlb to false.
}
\value{
These functions will invisibly return a three-column matrix, where the first two columns are interval start and end times and the third column is the number of taxa/lineages counted in that interval.
}
\author{David W. Bapst}
\seealso{
\code{
\link{multiDiv}},\code{\link{timeSliceTree}},

There are several different functions for traditional LTT plots (phylogenetic diversity curves), such as the function ,\code{\link{ltt.plot}} in the package ape, the function ltt() in the package phytools, the function plotLtt() in the package laser and the function LTT.average.root() in the package TreeSim.
}
\examples{
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's see what the 'true' diversity curve looks like in this case
#plot the FADs and LADs with taxicDivCont()
taxicDivCont(taxa)
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
  #using a tree from rtree in ape
set.seed(444)
tree<-rtree(100)
phyloDiv(tree)

#a neat example of using phyDiv with timeSliceTree 
 #to simulate doing molecular-phylogeny studies 
 #of diverification...in the past
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
taxicDivCont(taxa)
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve 
 #we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=TRUE)

#an example of a 'spiky' diversity curve and why split.int is a good thing
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
taxaDiv<-taxicDivCont(taxa)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
rangesDisc<-binTimeData(rangesCont,int.length=10)
#now let's plot with taxicDivDisc() but with the intervals from taxaDiv
	#by default, split.int=TRUE
taxicDivDisc(rangesDisc,int.times=taxaDiv[,1:2],split.int=TRUE)
#look pretty
#now let's turn off split.int
taxicDivDisc(rangesDisc,int.times=taxaDiv[,1:2],split.int=FALSE)
#looks 'spiky'!
}