\name{perCapitaRates}
\alias{perCapitaRates}
\title{perCapitaRates}
\description{Calculates and plots per-capita origination and extinction rates from sequential discrete-time taxon ranges, following Foote (2000).}
\usage{

perCapitaRates(


timeList, int.times = NULL, plot = TRUE, plotLogRich = FALSE, 
    timelims = NULL, extant.adjust = 0.01, split.int = TRUE)



}
\arguments{
  \item{timeList}{A list giving interval data. See details below.}


  \item{plot}{If true, diversity curve is plotted}
  \item{plotLogRich}{If true, taxic diversity plotted on log scale}

}
\details{
Background: Diversity curves are plots of species/taxon/lineage richness over time for a particular group of organisms. For paleontological studies, these are generally based on per-taxon range data while more recently in evolutionary biology, molecular phylogenies have been used to calculate lineage-through-time plots (LTTs). Neither of these approaches are without their particular weaknesses; reconstructing the true history of biodiversity is a difficult task no matter what data is available.

The diversity curves produced by these functions will always measure diversity within binned time intervals (and plot them as rectangular bins). For continuous-time data or phylogenies, one could decrease the int.length used to get what is essentially an 'instantaneous' estimate of diversity. This is warned against, however, as most historical diversity data will have some time-averaging or uncertain temporal resolution and thus is probably not finely-resolved enough to calculate instantaneous estimates of diversity.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

As diversity is counted within binned intervals, the true standing diversity may be somewhat lower than the measured (observed) quantity, particularly if intervals are longer than the mean duration of taxa is used. This will be an issue with all diversity curve functions, but particularly the discrete-time variant. For diversity data in particularly large discrete time intervals, plotting this data in smaller bins which do not line up completely with the original intervals will create a 'spiky' diversity curve, as these smaller intersecting bins will have a large number of taxa which may have been present in either of the neighboring intervals. This will give these small bins an apparently high estimated standing diversity. This artifact is avoided with the default setting split.int=TRUE, which will split any input or calculated intervals so that they start and end at the boundaries of the discrete-time range bins.

The timeList object should be a list composed of two matrices, the first matrix giving by-interval start and end times (in absolute time), the secont matrix giving the by-taxon first and last appearances in the intervals defined in the first matrix, numbered as the rows. Absolute time should be decreasing, while the intervals should be numbered so that the number increases with time. Unlike getSampProbDisc, the intervals can be overlapping. See the documentation for binTimeData for more information.

phyloDiv will resolve polytomies to be dichotomous nodes seperated by zero-length branches prior to calculating the diversity curve. There is no option to alter this behavior, but it should not affect the use of the function because the addition of the zero-length branches should produce an identical diveristy history as a polytomy. phyloDiv will also drop zero-length terminal branches, as with the function dropZLB. This the default behavior for the function but can be turned off by setting the argument drop.zlb to FALSE.
}
\value{
These functions will invisibly return a three-column matrix, where the first two columns are interval start and end times and the third column is the number of taxa/lineages counted in that interval.
}
\author{David W. Bapst}
\seealso{
\code{\link{DiversityCurves}}, \code{\link{binTimeData}}

There are several different functions for traditional LTT plots (phylogenetic diversity curves), such as the function ,\code{\link{ltt.plot}} in the package ape, the function \code{ltt} in the package phytools, the function \code{plotLtt} in the package laser and the function \code{LTT.average.root} in the package TreeSim.
}
\examples{
#Simulate some fossil ranges with simFossilTaxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=80,maxtaxa=100,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont <- sampleRanges(taxa,r=0.5)
#Now let's use binTimeData() to bin in intervals of 5 time units
rangesDisc <- binTimeData(rangesCont,int.length=5)
#and get the per-capita rates
perCapitaRates(rangesDisc)
#on a log scale
perCapitaRates(rangesDisc,plotLogRates=TRUE)

#get mean and median per-capita rates
res<-perCapitaRates(rangesDisc,plot=FALSE)
apply(res[,c("pRate","qRate")],2,mean,na.rm=TRUE)
apply(res[,c("pRate","qRate")],2,median,na.rm=TRUE)

#with modern taxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=80,maxtaxa=100,maxtime=1000)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont <- sampleRanges(taxa,r=0.5,,modern.samp.prob=1)
#Now let's use binTimeData() to bin in intervals of 5 time units
rangesDisc <- binTimeData(rangesCont,int.length=5)
#and now get per-capita rates
perCapitaRates(rangesDisc)
}