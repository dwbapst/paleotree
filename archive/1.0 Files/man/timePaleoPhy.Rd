\name{timePaleoPhy}
\alias{timePaleoPhy}
\alias{bin_timePaleoPhy}
\title{Timescaling of Paleo-Phylogenies}
\description{Takes an unscaled cladogram of fossil taxa and information on their ranges and outputs a timescaled tree using various methods. Also can resolve polytomies randomly and output samples or randomly-resolved trees.}
\usage{
timePaleoPhy(tree, timeData, type = "basic", vartime = NULL, ntrees = 1, randres = F, add.term = F, rand.obs = F, node.mins = NULL, plot = F)
bin_timePaleoPhy(tree, timeList, type = "basic", vartime = NULL, ntrees = 1, randres = F, sites = NULL,  add.term = F, rand.obs = F, node.mins = NULL, plot = F)
}
\arguments{
  \item{tree}{An unscaled cladogram of fossil taxa}
  \item{timeData}{Two-column matrix of first and last occurrances in absolute continous time, with rownames as the taxon IDs used on the tree}
  \item{type}{Type of time-scaling method used. Can be "basic", "equal", "aba", "zbla" or "mbl". See below for details. Type "basic" by default.}
  \item{vartime}{Time variable; usage depends on the method 'type' argument. Ignored if type = "basic"}
  \item{ntrees}{Number of time-scaled trees to output}
  \item{randres}{Should polytomies be randomly resolved?}
  \item{add.term}{If true, adds terminal ranges (see below)}
  \item{rand.obs}{Should the tips represent observation times uniform distributed within taxon ranges?}
  \item{node.mins}{Minimum ages of nodes on the tree, see below}
  \item{plot}{If true, plots the input and output phylogenies}
  \item{timeList}{A list composed of two matrices giving interval times and taxon appearance datums, as would be output by binTimeData. The rownames of the second matrix should be the taxon IDs}
  \item{sites}{A two column matrix, composed of site IDs for taxon FADs and LADs. Does not have to be given by default; see explanation below.}
}
\details{
This function is an attempt to unify and collect previously used and discussed methods for time-scaling phylogenies of fossil taxa. Unfortunately, it is difficult to attribute some method types to specific references.

There are five method types that can be used by timePaleoPhy. Four of these use some amount of absolute time, chosen a priori, to time-scale the tree. This is handled by the argument vartime, which is NULL by default.

	"basic" : This most simple of methods ignores vartime and scales nodes so they are as old as the first appearance of their oldest descendant. This method produces many zero-length branches.

	"equal" : The 'equal' method defined by G. Lloyd and used in Brusatte et al. (2008) and Lloyd et al. (2012). Originally usable in code supplied by G. Lloyd, it is recreated here. This method works by increasing the time of the root divergence by some amount (here set by vartime) and then adjusting zero-length branches so that time on early branches is re-apportioned out along those later branches equally.

	"aba" : All branches additive. This method takes the "basic" tree and adds vartime to all branches.

	"zlba" : Zero-length branches additive. This method adds vartime to all zero-length branches in the "basic" tree. Mentioned by Hunt and Carrano, 2010.

	"mbl": Minimum branch length. Scales all branches so they are greater than or equal to vartime, and subtract time added to later branches from earlier branches in order to maintain the temporal structure of events.

This function cannot time-scale branches relative to reconstructed character changes along branches, as used by Brusatte et al. (2008).

By default, timePaleoPhy does not resolve polytomies, instead outputting a time-scaled tree that is only as resolved as the input tree. If randres=T, then polytomies will be randomly resolved using multi2di() from the package ape. If randres=T and ntrees=1, then a warning is printed that users should analyze multiple randomly-resolved trees, rather than a single such tree. If ntrees is greater than one and randres=F, the function will fail and a warning is issued, as these are by definition arguments in conflict with each other.

By default, this function will not add the ranges of taxa when time-scaling a tree, so that the tips correspond temporally to the first appearance datums of the given taxa. If add.term=T, then the 'terminal ranges' of the taxa are added to the tips after tree is time-scaled, such that the tips now correspond to the last appearance datums. 

If add.term and rand.obs are both set to be true, then it is assumed that users wish the tips to represent observations made with some temporal uncertainty, such that they might have come from any point within a taxon's range. This might be the case, for example, if a user is interested in applying phylogeny-based approaches to studying trait evolution, but have per-taxon measurements of traits that come from museum specimens with uncertain temporal placement. When these arguments are set to be true, the tips are placed randomly within taxon ranges, as if uniformly distributed. Obviously, as with randres, multiple trees should be created and then analyzed. If add.term=F and rand.obs=T, the function fails and a warning is issued.

The minimum dates of nodes can be set using node.mins; this argument takes a vector of the same length as the number of nodes, with dates given in the same order as nodes are number in the tree$edge matrix. Not all nodes need be set; those without minimum dates can be given as NA in node.mins.

All trees are output with an element $root.time. This is the time of the root on the tree and is important for comparing patterns across trees.

These functions will intuitively drop taxa from the tree with NA for range or are missing from timeData or timeList.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

bin_timePaleoPhy is a wrapper of timePaleoPhy which produces timescaled trees for datasets which only have interval data available. For each output tree, taxon FADs and LADs are placed within their listed intervals under a uniform distribution. Thus, a large sample of time-scaled trees will approximate the uncertainty in the actual timing of the FADs and LADs.

The sites argument allows users to constrain the placement of dates in bin_timePaleoPhy by restricting multiple fossil taxa whose FADs or LADs are from the same very temporally restricted sites (such as fossil-rich Lagerstatten) to always have the same date, across many iterations of time-scaled trees from bin_timePaleoPhy. To do this, simply give a matrix where the "site" of each FAD and LAD for every taxon is listed, as corresponding to the second matrix in timeList. If no sites matrix is given (the default), then it is assumed all fossil come from different "sites" and there is no shared temporal structure among the events.
}
\value{
The output of these functions is a time-scaled tree or set of time-scaled trees, of either class phylo or multiphylo, depending on the argument ntrees.
}
\references{
Brusatte, S. L., M. J. Benton, M. Ruta, and G. T. Lloyd. 2008. Superiority, Competition, and Opportunism in the Evolutionary Radiation of Dinosaurs. Science 321(5895):1485-1488.

Hunt, G., and M. T. Carrano. 2010. Models and methods for analyzing phenotypic evolution in lineages and clades. In J. Alroy, and G. Hunt, eds. Short Course on Quantitative Methods in Paleobiology. Paleontological Society.

Lloyd, G. T., S. C. Wang, and S. L. Brusatte. 2012. Identifying Heterogeneity in Rates of Morphological Evolutio: Discrete Character Change in the Evolution of Lungfish(Sarcopterygii, Dipnoi). Evolution 66(2):330-348.
}
\author{David W. Bapst, heavily inspired by code supplied by Graeme Lloyd and Gene Hunt.}

\seealso{
\code{\link{srcTimePaleoPhy}},\code{\link{binTimeData}},\code{\link{multi2di}}
}
\examples{
#timePaleoPhy
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy() using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)

#that tree didn't look very resolved, does it? (The curse of sampled ancestry!)
#if we set ntrees>1, timePaleoPhy() will make multiple time-trees
#each resulting tree will have polytomies randomly resolved in different ways using multi2di()
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=1,randres=TRUE,add.term=TRUE,plot=TRUE)
#notice that the warning it prints!
#now let's plot the first tree (both trees will be identical because we used set.seed)
phyloDiv(ttree)
#we would need to set ntrees to a large number to get a fair sample of trees

#compare different methods of timePaleoPhy
layout(matrix(1:6,3,2));par(mar=c(3,2,1,2))
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="basic",vartime=NULL,add.term=TRUE)))
	axisPhylo();text(x=50,y=23,"type=basic",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="equal",vartime=10,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=equal",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="aba",vartime=1,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=aba",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="zlba",vartime=1,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=zlba",adj=c(0,0.5),cex=1.2)
plot(ladderize(timePaleoPhy(cladogram,rangesCont,type="mbl",vartime=1,add.term=TRUE)));axisPhylo()
	axisPhylo();text(x=55,y=23,"type=mbl",adj=c(0,0.5),cex=1.2)

#Using bin_timePaleoPhy to timescale with discrete interval data
#first let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
ttree<-bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,randres=TRUE,add.term=TRUE,plot=TRUE)
#notice that the warning it prints!
phyloDiv(ttree)
}
\keyword{datagen}
