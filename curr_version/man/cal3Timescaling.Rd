\name{cal3Timescaling}
\alias{cal3TimePaleoPhy}
\alias{bin_cal3TimePaleoPhy}
\title{Three Rate Calibrated Timescaling of Paleo-Phylogenies}
\description{
Timescales an unscaled cladogram of fossil taxa, using information on their ranges and estimates of the instantaneous rates of branching, extinction and sampling. The output is a sample of timescaled trees, as resulting from a stochastic algorithm which samples observed gaps in the fossil record with weights calculated based on the input rate estimates. This function also uses the three-rate calibrated time-scaling algorithm to stochastically resolve polytomies and infer potential ancestor-descendant relationships, simultaneous with the time-scaling treatment.}
\usage{
cal3TimePaleoPhy(tree, timeData, brRate, extRate, sampRate, ntrees = 1, anc.wt = 1,  
    node.mins = NULL, rand.obs = FALSE, FAD.only = FALSE, adj.obs.wt = TRUE,
    root.max = 200, step.size = 0.1, plot = FALSE)

bin_cal3TimePaleoPhy(tree, timeList, brRate, extRate, sampRate, ntrees = 1, 
    nonstoch.bin=FALSE,sites = NULL, anc.wt = 1, node.mins = NULL, rand.obs = FALSE, 
    FAD.only = FALSE, adj.obs.wt = TRUE, root.max = 200, step.size = 0.1, plot = FALSE)
}
\arguments{
  \item{tree}{An unscaled cladogram of fossil taxa}

  \item{timeData}{Two-column matrix of first and last occurrances in absolute continous time, with rownames as the taxon IDs used on the tree}

  \item{brRate}{Either a single estimate of the instanteous rate of branching (also known as the 'per-capita' origination rate, or speciation rate if taxonomic level of interest is species) or a vector of per-taxon estimates}

  \item{extRate}{Either a single estimate of the instanteous extinction rate (also known as the 'per-capita' extinction rate) or a vector of per-taxon estimates}

  \item{sampRate}{Either a single estimate of the instanteous sampling rate or a vector of per-taxon estimates}

  \item{ntrees}{Number of time-scaled trees to output}

  \item{anc.wt}{Weighting against inferring ancestor-descendant relationships. The argument anc.wt allows users to change the default consideration of anc-desc relationships. This value is used as a multiplier applied to the probability of choosing any node position which would infer an ancestor-descendant relationship. By default, anc.wt=1, and thus these probabilities are unaltered. if anc.wt is less than 1, the probabilities decrease and at anc.wt=0, no ancestor-descendant relationships are inferred at all. Can be a single value or a vector of per-taxon values, such as if a user wants to only allow plesiomorphic taxa to be ancestors.}

  \item{rand.obs}{Should the tips represent observation times uniform distributed within taxon ranges? If rand.obs is TRUE, then it is assumed that users wish the tips to represent observations made with some temporal uncertainty, such that they might have come from any point within a taxon's range. This might be the case, for example, if a user is interested in applying phylogeny-based approaches to studying trait evolution, but have per-taxon measurements of traits that come from museum specimens with uncertain temporal placement. When rand.obs=T, the tips are placed randomly within taxon ranges, as if uniformly distributed.}

  \item{node.mins}{Minimum ages of nodes on the tree. The minimum dates of nodes can be set using node.mins; this argument takes a vector of the same length as the number of nodes, with dates given in the same order as nodes are they are numbered in the tree$edge matrix (note that in tree$edge, the tips are given first Ntip numbers and these are ignored here). Not all nodes need be set; those without minimum dates can be given as NA in node.mins. These nodes will be frozen and will not be shifted by the cal3 algorithm. If the dates refer to a polytomy, then the first divergence will be frozen with additional divergence able to occur after the minimum date.}

  \item{FAD.only}{Should the tips represent observation times at the start of the taxon ranges? If TRUE, result is similar to when terminal ranges are no added on with timePaleoPhy. If FAD.only and rand.obs are both TRUE, the function will stop and a warning will be produced.}

  \item{adj.obs.wt}{If the time of observation are before the LAD of a taxon, should the weight of the time of observation be adjusted to account for the known observed history of the taxon which occurs AFTER the time of observation? Should only have an effect if time of observation ISN'T the LAD, if times of observation for a potential ancestor are earlier than the first apperance of their potential descendents and if the ancestral weights for taxa aren't set to zero (so there can be potential ancestors).}

  \item{root.max}{Maximum time before the first FAD that the root can be pushed back to}

  \item{step.size}{Step size of increments used in zipper algorithm to assign node ages}

  \item{plot}{If true, plots the input, "basic" timescaled and output cal3-timescaled phylogenies}

  \item{timeList}{A list composed of two matrices giving interval times and taxon appearance datums, as would be output by binTimeData. The rownames of the second matrix should be the taxon IDs}

  \item{nonstoch.bin}{If true, dates are not stochastically pulled from uniform distributions. See below for more details.}

  \item{sites}{Optional two column matrix, composed of site IDs for taxon FADs and LADs. The sites argument allows users to constrain the placement of dates by restricting multiple fossil taxa whose FADs or LADs are from the same very temporally restricted sites (such as fossil-rich Lagerstatten) to always have the same date, across many iterations of time-scaled trees. To do this, simply give a matrix where the "site" of each FAD and LAD for every taxon is listed, as corresponding to the second matrix in timeList. If no sites matrix is given (the default), then it is assumed all fossil come from different "sites" and there is no shared temporal structure among the events.}
}
\details{
The three-rate calibrated ("cal3") algorithm time-scales trees by stochastically picking node divergence times relative to a probability distribution of expected waiting times between speciation and first appearance in the fossil record. This algorithm is extended to apply to resolving polytomies and designating possible ancestor-descendant relationships. The full details of this method, its sister method the sampling-rate-calibrated timescaling method and the details of the algorithm will be given in Bapst (in prep for Paleobiology).

Briefly, cal3 timescaling is done by examining each node seperately, moving from the root upwards. Ages of branching nodes are constained below by the ages of the nodes below them (except the root; hence the need for the root.max argument) and constrained above by the first appearance dates (FADs) of the daughter lineages. The position of the branching event within this constained range implies different amounts of unobserved evolutionary history. cal3 considers a large number of potential positions for the branching node (the notation in the code uses the analogy of viewing the branching event as a 'zipper') and calculates the summed unobservered evolutionary history implied by each branching time. The probability density of each position is then calculated under a gamma distribution with a shape parameter of 2 (implying that it is roughly the sum of two normal waiting times under an exponential) and a rate parameter which takes into account both the probability of not observing a lineage of a certain duration and the 'twiginess' of the branch, i.e. the probability of having short-lived descendants which went extinct and never were sampled (ala Friedman and Brazeau, 2011). These densities calculated under the gamma distribution are then used as weights to stochastically sample the possible positions for the branching node. This basic framework is extended to polytomies by allowing a branching event to fall across multiple potential lineages, adding each lineage one by one, from earliest appearing to latest appearing (the code notation refers to this as a 'parallel zipper').

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

These functions will intuitively drop taxa from the tree with NA for range or that are missing from timeData.

The sampling rate used by cal3 methods is the instantaneous sampling rate, as estimated by various other function in the paleotree package. See getSampRateCont for more details. If you have the per-time unit sampling probability ('R' as opposed to 'r') look at the sampling parameter conversion functions also included in this package. Most datasets will probably use getSampProbDisc and sProb2sRate prior to using this function, as shown in an example below.

The branching and extinction rate are the 'per-capita' instantaneous origination/extinction rates for the taxic level of the tips of the tree being time-scaled. Any user of the cal3 timescaling method has multiple options for estimating these rates. One is to seperately calculate the per-capita rates (following the equations in Foote, 2001) across multiple intervals and take the mean for each rate. A second, less preferred option, would be to use the extinction rate calculated from the sampling rate above (under ideal conditions, this should be very close to the mean 'per-capita' rate calculated from by-interval FADs and LADs). The branching rate in this case could be assumed to be very close to the extinction rate, given the tight relationship observed in general between these two (Stanley, 1976; see Foote et al., 1999, for a defense of this approach), and thus the extinction rate estimate could be used also for the branching rate estimate. (This is what is done for the examples below.) A third option for calculating all three rates simultaneously would be to apply likelihood methods developed by Foote (2002) to forward and reverse survivorship curves. Note that only one of these three suggested methods is implemented in paleotree: estimating the sampling and extinction rates from the distribution of taxon durations via getSampRateCont and getSampProbDisc.

By default, the cal3 functions will consider that ancestor-descendant relationships may exist among the given taxa, under a budding cladogenetic or anagenetic modes. Which tips are designated as which is given by two additional elements added to the output tree, $budd.tips (taxa designated as ancestors via budding cladogenesis) and $anag.tips (taxa designated as ancestors via anagenesis). This can be turned off by setting anc.wt=0. As this function may infer anagenetic relationships during time-scaling, this can create zero-length terminal branches in the output. Use \code{\link{dropZLB}} to get rid of these before doing analyses of lineage diversification.

Unlike timePaleoPhy, cal3 methods will always resolve polytomies. In general, this is done using the rate calibrated algorithm, although if argument randres=TRUE, polytomies will be randomly resolved with uniform probability, ala multi2di from ape. Also, cal3 will always add the terminal ranges of taxa. However, because of the ability to infer potential ancestor-descendant relationships, the length of terminal branches may be shorter than taxon ranges themselves, as budding may have occurred during the range of a morphologically static taxon. By resolving polytomies with the cal3 method, this function allows for taxa to be ancestral to more than one descendant taxon. Thus, users who believe their dataset may contain indirect ancestors are encouraged by the package author to try cal3 methods with their consensus trees, as opposed to using the set of most parsimonious trees. Comparing the results of these two approaches may be very revealing.

cal3TimePaleoPhy is only applicable to datasets with taxon occurances in continuous time. In comparison, bin_cal3TimePaleoPhy is a wrapper of cal3TimePaleoPhy which produces timescaled trees for datasets which only have interval data available. For each output tree, taxon FADs and LADs are placed within their listed intervals under a uniform distribution. Thus, a large sample of time-scaled trees will approximate the uncertainty in the actual timing of the FADs and LADs.

By setting the argument nonstoch.bin to TRUE in bin_cal3TimePaleoPhy, the dates are NOT stochastically pulled from uniform bins but instead FADs are assigned to the earliest time of whichever interval they were placed in and LADs are placed at the most recent time in their placed interval. This option may be useful for plotting. The sites argument becomes arbitrary if nonstoch.bin is TRUE.

If timeData or the elements of timeList are actually data.frames (as output by read.csv or read.table), these will be coerced to a matrix.
}
\note{
Most importantly, please note the stochastic element of the three rate-calibrated time-scaling methods. These do not use traditional optimization methods, but instead draw divergence times from a distribution defined by the probability of intervals of unobserved evolutionary history. This means analyses MUST be done over many cal3-timescaled trees for analytical rigor! No one tree is correct.
}
\value{
The output of these functions is a time-scaled tree or set of time-scaled trees, of either class phylo or multiphylo, depending on the argument ntrees. All trees are output with an element $root.time. This is the time of the root on the tree and is important for comparing patterns across trees.

Trees created with bin_cal3TimePaleoPhy will output with some additional elements, in particular $ranges.used, a matrix which records the continuous-time ranges generated for time-scaling each tree. (Essentially a pseudo-timeData matrix.)
}
\references{
Bapst, in prep. Time-scaling Trees of Fossil Taxa. To be submitted to \emph{Paleobiology}

Foote, M. 2000. Origination and extinction components of taxonomic diversity: general problems. Pp. 74-102. In D. H. Erwin, and S. L. Wing, eds. Deep Time: Paleobiology's Perspective. The Paleontological Society, Lawrence, Kansas.

Foote, M. 2001. Inferring temporal patterns of preservation, origination, and extinction from taxonomic survivorship analysis. \emph{Paleobiology} 27(4):602-630.

Friedman, M., and M. D. Brazeau. 2011. Sequences, stratigraphy and scenarios: what can we say about the fossil record of the earliest tetrapods? \emph{Proceedings of the Royal Society B: Biological Sciences} 278(1704):432-439.

Stanley, S. M. 1979. Macroevolution: Patterns and Process. W. H. Freeman, Co., San Francisco.
}
\author{David W. Bapst}
\seealso{
\code{\link{timePaleoPhy}}, \code{\link{binTimeData}}, \code{\link{getSampRateCont}}, \code{\link{multi2di}}
}
\examples{
#Simulate some fossil ranges with simFossilTaxa
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges
rangesCont <- sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram <- taxa2cladogram(taxa,plot=TRUE)
#this library allows one to use rate calibrated type time-scaling methods (Bapst, in prep.)
#to use these, we need an estimate of the sampling rate (we set it to 0.5 above)
SRres <- getSampRateCont(rangesCont)
sRate <- SRres[[2]][2]
#we also need extinction rate and branching rate; we can get extRate from getSampRateCont too
#we'll assume extRate=brRate (ala Foote et al., 1999); may not always be a good assumption
divRate<-SRres[[2]][1]
#now let's try cal3TimePaleoPhy, which timescales using a sampling rate to calibrate
#This can also resolve polytomies based on sampling rates, with some stochastic decisions
ttree <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    sampRate=sRate,ntrees=1,plot=TRUE)
#notice the warning it gives!
phyloDiv(ttree)

#by default, cal3TimePaleoPhy is allowed to predict indirect ancestor-descendant relationships
#can turn this off by setting anc.wt=0
ttree <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    sampRate=sRate,ntrees=1,anc.wt=0,plot=TRUE)

#let's look at how three trees generated with very different time of obs. look
ttreeFAD <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    FAD.only=TRUE,rand.obs=FALSE,sampRate=sRate,ntrees=1,plot=TRUE)
ttreeRand <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    FAD.only=FALSE,rand.obs=TRUE,sampRate=sRate,ntrees=1,plot=TRUE)
#by default the time of observations are the LADs
ttreeLAD <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    FAD.only=FALSE,rand.obs=FALSE,sampRate=sRate,ntrees=1,plot=TRUE)
layout(1:3);parOrig<-par(mar=c(0,0,0,0))
plot(ladderize(ttreeFAD));plot(ladderize(ttreeRand));plot(ladderize(ttreeLAD))

#to get a fair sample of trees, let's increase ntrees
ttrees <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    sampRate=sRate,ntrees=9,plot=FALSE)
#let's compare nine of them at once in a plot
layout(matrix(1:9,3,3));parOrig<-par(mar=c(0,0,0,0))
for(i in 1:9){plot(ladderize(ttrees[[i]]),show.tip.label=FALSE)}
#they are all a bit different!

#can plot the median diversity curve with multiDiv
layout(1); par(parOrig)
multiDiv(ttrees)

#using node.mins
#let's say we have (molecular??) evidence that node #5 is at least 1200 time-units ago
nodeDates <- rep(NA,(Nnode(cladogram)-1))
nodeDates[5]<-1200
ttree <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
    sampRate=sRate,ntrees=1,node.mins=nodeDates,plot=TRUE)

#example with time in discrete intervals
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont <- sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram <- taxa2cladogram(taxa,plot=TRUE)
#Now let's use binTimeData to bin in intervals of 1 time unit
rangesDisc <- binTimeData(rangesCont,int.length=1)
#we can do something very similar for the discrete time data (can be a bit slow)
SPres <- getSampProbDisc(rangesDisc)
sProb <- SPres[[2]][2]
#but that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
sRate1 <- sProb2sRate(sProb,int.length=1)
#we also need extinction rate and branching rate (see above)
divRate<-SPres[[2]][1]
#estimates that r=0.3... kind of low (simulated sampling rate is 0.5)
#Note: for real data, you may need to use an average int.length (no constant length)
ttree <- bin_cal3TimePaleoPhy(cladogram,rangesDisc,brRate=divRate,extRate=divRate,
    sampRate=sRate1,ntrees=1,plot=TRUE)
phyloDiv(ttree)
#can also force the appearance timings not to be chosen stochastically
ttree1 <- bin_cal3TimePaleoPhy(cladogram,rangesDisc,brRate=divRate,extRate=divRate,
    sampRate=sRate1,ntrees=1,nonstoch.bin=TRUE,plot=TRUE)
phyloDiv(ttree1)


#example with multiple values of anc.wt
ancWt <- sample(0:1,nrow(rangesDisc[[2]]),replace=TRUE)
names(ancWt)<-rownames(rangesDisc[[2]])
ttree1 <- bin_cal3TimePaleoPhy(cladogram,rangesDisc,brRate=divRate,extRate=divRate,
    sampRate=sRate1,ntrees=1,anc.wt=ancWt,plot=TRUE)
}