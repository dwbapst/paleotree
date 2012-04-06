\name{srcTimescaling}
\alias{srcTimePaleoPhy}
\alias{bin_srcTimePaleoPhy}
\title{SampRate-Calibrated Timescaling of Paleo-Phylogenies}
\description{
This function takes an input an unscaled cladogram of fossil taxa, information on their ranges and an estimate of the instantaneous rate of sampling. The output is a sample of timescaled trees, as resulting from a stochastic algorithm that samples observed gaps in the fossil record with weights calculated from the sampling rate. This function also uses the sampling-rate calibrated time-scaling algorithim to resolve polytomies randomly and infer potential ancestor-descendant relationships, simultaneous with the time-scaling.}
\usage{
srcTimePaleoPhy(tree, timeData, sampRate, ntrees = 1, anc.wt = 1, rand.obs = F, node.mins = NULL,
root.max = 200, plot = F)

bin_srcTimePaleoPhy(tree, timeList, sampRate, ntrees = 1, sites = NULL, anc.wt = 1, node.mins = NULL,
 rand.obs = F, root.max = 200, plot = F)
}
\arguments{
  \item{tree}{An unscaled cladogram of fossil taxa}

  \item{timeData}{Two-column matrix of first and last occurrances in absolute continous time, with rownames as the taxon IDs used on the tree}

  \item{sampRate}{Either a single estimate of the instanteous sampling rate or a vector of per-taxon estimates}

  \item{ntrees}{Number of time-scaled trees to output}

  \item{anc.wt}{Weighting against inferring ancestor-descendant relationships. The argument anc.wt allows users to change the default consideration of anc-desc relationships. This value is used as a multiplier applied to the probability of choosing any node position which would infer an ancestor-descendant relationship. By default, anc.wt=1, and thus these probabilities are unaltered. if anc.wt is less than 1, the probabilities decrease and at anc.wt=0, no ancestor-descendant relationships are inferred at all.}

  \item{rand.obs}{Should the tips represent observation times uniform distributed within taxon ranges? If rand.obs=T, then it is assumed that users wish the tips to represent observations made with some temporal uncertainty, such that they might have come from any point within a taxon's range. This might be the case, for example, if a user is interested in applying phylogeny-based approaches to studying trait evolution, but have per-taxon measurements of traits that come from museum specimens with uncertain temporal placement. When rand.obs=T, the tips are placed randomly within taxon ranges, as if uniformly distributed.}

  \item{node.mins}{Minimum ages of nodes on the tree. The minimum dates of nodes can be set using node.mins; this argument takes a vector of the same length as the number of nodes, with dates given in the same order as nodes are number in the tree$edge matrix. Not all nodes need be set; those without minimum dates can be given as NA in node.mins. These nodes will be frozen and will not be shifted by the SRC algorithm. If the dates refer to a polytomy, then the first divergence will be frozen with additional divergence able to occur after the minimum date.}

  \item{root.max}{Maximum time before the first FAD that the root can be pushed back to}

  \item{plot}{If true, plots the input, "basic" timescaled and output SRC-timescaled phylogenies}

  \item{timeList}{A list composed of two matrices giving interval times and taxon appearance datums, as would be output by binTimeData. The rownames of the second matrix should be the taxon IDs}

  \item{sites}{Optional two column matrix, composed of site IDs for taxon FADs and LADs. The sites argument allows users to constrain the placement of dates in bin_srcTimePaleoPhy by restricting multiple fossil taxa whose FADs or LADs are from the same very temporally restricted sites (such as fossil-rich Lagerstatten) to always have the same date, across many iterations of time-scaled trees from bin_timePaleoPhy. To do this, simply give a matrix where the "site" of each FAD and LAD for every taxon is listed, as corresponding to the second matrix in timeList. If no sites matrix is given (the default), then it is assumed all fossil come from different "sites" and there is no shared temporal structure among the events.}
}
\details{
The sampling-rate calibrated (SRC) algorithim time-scales trees by stochastically picking node divergence times relative to a probability distribution of expected waiting times between speciation and first appearance in the fossil record. This algorithm is also extended to apply to resolving polytomies and designating possible ancestor-descendant relationships. The full details of this method and the algorithm use will be given in Bapst (in prep). Its performance with other time-scaling methods will also be compared via simulation.

srcTimePaleoPhy is only applicable to datasets with taxon occurances in continuous time. bin_srcTimePaleoPhy is a wrapper of srcTimePaleoPhy which produces timescaled trees for datasets which only have interval data available. For each output tree, taxon FADs and LADs are placed within their listed intervals under a uniform distribution. Thus, a large sample of time-scaled trees will approximate the uncertainty in the actual timing of the FADs and LADs.

The sampling rate used by SRC methods is the instantaneous sampling rate, as estimated by various other function in the paleotree package. See getSampRateCont for more details. If you have the per-time unit sampling probability ('R' as opposed to 'r') look at the sampling parameter conversion functions also included in this package. Most datasets will probably use getSampProbDisc and sProb2sRate prior to using this function, as shown in an example below.

By default, the SRC functions will consider that ancestor-descendant relationships may exist among the given taxa, under a budding cladogenetic or anagenetic modes. Which tips are designated as which is given by two additional elements added to the output tree, $budd.tips (taxa designated as ancestors via budding cladogenesis) and $anag.tips (taxa designated as ancestors via anagenesis). This can be turned off by setting anc.wt=0. As this function may infer anagenetic relationships during time-scaling, this can create zero-length terminal branches in the output. Use dropZLB() to get rid of these before doing analyses of lineage diversification.

Unlike timePaleoPhy, SRC methods will always resolve polytomies (using the sampling-rate calibrated algorithim) and will always add the terminal ranges of taxa. However, because of the ability to infer potential ancestor-descendant relationships, the length of terminal branches may be shorter than taxon ranges themselves, as budding may have occurred during the range of a morphologically static taxon. By resolving polytomies with the SRC method, this function allows for taxa to be ancestral to more than one descendant taxon.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

These functions will intuitively drop taxa from the tree with NA for range or that are missing from timeData.
}
\note{
Most importantly, please note the stochastic element of the SRC method. It does not use traditional optimization methods, but instead pulls node times from a distribution. This means analyses MUST be done over many SRC-timescaled trees for analytical rigor! No one tree is correct.
}
\value{
The output of these functions is a time-scaled tree or set of time-scaled trees, of either class phylo or multiphylo, depending on the argument ntrees. All trees are output with an element $root.time. This is the time of the root on the tree and is important for comparing patterns across trees.
}
\references{
Bapst, in prep. Time-scaling Trees of Fossil Taxa. To be submitted to Paleobiology.
}
\author{David W. Bapst}
\seealso{
\code{\link{timePaleoPhy}},\code{\link{binTimeData}},\code{\link{getSampRateCont}},\code{\link{multi2di}}
}
\examples{
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#this library allows one to use SRC type time-scaling methods (Bapst, in prep.)
#to use these, we need an estimate of the sampling rate (we set it to 0.5 above)
SRres<-getSampRateCont(rangesCont)
sRate<-SRres$pars[2]
#now let's try srcTimePaleoPhy(), which timescales using a sampling rate to calibrate
#This can also resolve polytomies based on sampling rates, with some stochastic decisions
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,plot=TRUE)
#notice the warning it gives!
phyloDiv(ttree)

#by default, srcTimePaleoPhy() is allowed to predict indirect ancestor-descendant relationships
#can turn this off by setting anc.wt=0
ttree<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=1,anc.wt=0,plot=TRUE)

#to get a fair sample of trees, let's increse ntrees
#can do an example of such an analysis via multDiv
ttrees<-srcTimePaleoPhy(cladogram,rangesCont,sampRate=sRate,ntrees=10,plot=FALSE)
multiDiv(ttrees)

#example with time in discrete intervals
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
#we can do something very similar for the discrete time data (can be a bit slow)
SPres<-getSampProbDisc(rangesDisc)
sProb<-SPres$pars[2]
#but that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
sRate1<-sProb2sRate(sProb,int.length=1)
#estimates that r=0.3... kind of low (simulated sampling rate is 0.5)
#Note: for real data, you may need to use an average int.length (no constant length)
ttree<-bin_srcTimePaleoPhy(cladogram,rangesDisc,sampRate=sRate1,ntrees=1,plot=TRUE)
phyloDiv(ttree)
}
\keyword{datagen}