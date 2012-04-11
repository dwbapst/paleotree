\name{timePaleoPhy}
\alias{timePaleoPhy}
\alias{bin_timePaleoPhy}
\title{Timescaling of Paleo-Phylogenies}
\description{Takes an unscaled cladogram of fossil taxa and information on their ranges and outputs a timescaled tree using various methods. Also can resolve polytomies randomly and output samples or randomly-resolved trees.}
\usage{
timePaleoPhy(tree, timeData, type = "basic", vartime = NULL, ntrees = 1, randres = FALSE, 
	add.term = FALSE, rand.obs = FALSE, node.mins = NULL, plot = FALSE)

bin_timePaleoPhy(tree, timeList, type = "basic", vartime = NULL, ntrees = 1,  nonstoch.bins=FALSE, 
	randres = FALSE, sites = NULL,  add.term = FALSE, rand.obs = FALSE, node.mins = NULL, plot = FALSE)
}
\arguments{
  \item{tree}{An unscaled cladogram of fossil taxa}

  \item{timeData}{Two-column matrix of first and last occurrances in absolute continous time, with rownames as the taxon IDs used on the tree}

  \item{type}{Type of time-scaling method used. Can be "basic", "equal", "aba", "zbla" or "mbl". Type="basic" by default. See the note below for more details}

  \item{vartime}{Time variable; usage depends on the method 'type' argument. Ignored if type = "basic"}

  \item{ntrees}{Number of time-scaled trees to output. If ntrees is greater than one and both randres and rand.obs are false, the function will fail and a warning is issued, as these arguments would simply produce multiple identical time-scaled trees.}

  \item{randres}{Should polytomies be randomly resolved? By default, timePaleoPhy does not resolve polytomies, instead outputting a time-scaled tree that is only as resolved as the input tree. If randres=T, then polytomies will be randomly resolved using multi2di() from the package ape. If randres=T and ntrees=1, a warning is printed that users should analyze multiple randomly-resolved trees, rather than a single such tree, although a tree is still output.}

  \item{add.term}{If true, adds terminal ranges.By default, this function will not add the ranges of taxa when time-scaling a tree, so that the tips correspond temporally to the first appearance datums of the given taxa. If add.term=T, then the 'terminal ranges' of the taxa are added to the tips after tree is time-scaled, such that the tips now correspond to the last appearance datums.}

  \item{rand.obs}{Should the tips represent observation times uniform distributed within taxon ranges? If rand.obs=T but add.term=F, the function fails and a warning is issued. If rand.obs=TRUE, tips are placed randomly within taxon ranges, as if uniformly distributed. This serves those users that wish for tips to represent observations made with some temporal uncertainty, such that they might have come from any point within a taxon's range. This might be the case, for example, if a user is interested in applying phylogeny-based approaches to studying trait evolution, but have per-taxon measurements of traits that come from museum specimens with uncertain temporal placement. As with randres, multiple trees should be created and then analyzed. }

  \item{node.mins}{Minimum ages of nodes on the tree. The minimum dates of nodes can be set using node.mins; this argument takes a vector of the same length as the number of nodes, with dates given in the same order as nodes are they are numbered in the tree$edge matrix (note that in tree$edge, the tips are given first Ntip numbers and these are ignored here). Not all nodes need be set; those without minimum dates can be given as NA in node.mins.}

  \item{plot}{If true, plots the input and output phylogenies}

  \item{timeList}{A list composed of two matrices giving interval times and taxon appearance datums, as would be output by binTimeData. The rownames of the second matrix should be the taxon IDs}

  \item{nonstoch.bin}{If true, dates are not stochastically pulled from uniform distributions. See below for more details.}

  \item{sites}{Optional two column matrix, composed of site IDs for taxon FADs and LADs. The sites argument allows users to constrain the placement of dates in bin_timePaleoPhy by restricting multiple fossil taxa whose FADs or LADs are from the same very temporally restricted sites (such as fossil-rich Lagerstatten) to always have the same date, across many iterations of time-scaled trees from bin_timePaleoPhy. To do this, simply give a matrix where the "site" of each FAD and LAD for every taxon is listed, as corresponding to the second matrix in timeList. If no sites matrix is given (the default), then it is assumed all fossil come from different "sites" and there is no shared temporal structure among the events.}
}
\details{
This function is an attempt to unify and collect previously used and discussed methods for time-scaling phylogenies of fossil taxa. Unfortunately, it is difficult to attribute some time-scaling methods to specific references in the literature. A review of these time-scaling methods is forthcoming (Bapst, in prep).

There are five method types that can be used by timePaleoPhy. Four of these use some value of absolute time, chosen a priori, to time-scale the tree. This is handled by the argument vartime, which is NULL by default and unused for type "basic".

\describe{
\item{"basic"}{This most simple of methods ignores vartime and scales nodes so they are as old as the first appearance of their oldest descendant. This method produces many zero-length branches (Smith, 1994).}
\item{"equal"}{The 'equal' method defined by G. Lloyd and used in Brusatte et al. (2008) and Lloyd et al. (2012). Originally usable in code supplied by G. Lloyd, it is recreated here. This method works by increasing the time of the root divergence by some amount (here set by vartime) and then adjusting zero-length branches so that time on early branches is re-apportioned out along those later branches equally.}
\item{"aba"}{All branches additive. This method takes the "basic" tree and adds vartime to all branches.}
\item{"zlba"}{Zero-length branches additive. This method adds vartime to all zero-length branches in the "basic" tree. Discussed by Hunt and Carrano, 2010.}
\item{"mbl"}{Minimum branch length. Scales all branches so they are greater than or equal to vartime, and subtract time added to later branches from earlier branches in order to maintain the temporal structure of events.}
}

These functions cannot time-scale branches relative to reconstructed character changes along branches, as used by Lloyd et al. (2012).

These functions will intuitively drop taxa from the tree with NA for range or are missing from timeData or timeList.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.

timePaleoPhy can only be applied to datasets where taxon appearances are in continuous time. bin_timePaleoPhy is a wrapper of timePaleoPhy which produces timescaled trees for datasets which only have interval data available. For each output tree, taxon FADs and LADs are placed within their listed intervals under a uniform distribution. Thus, a large sample of time-scaled trees will approximate the uncertainty in the actual timing of the FADs and LADs. 

By setting the argument nonstoch.bin to TRUE, the dates are NOT stochastically pulled from uniform bins but instead FADs are assigned to the earliest time of whichever interval they were placed in and LADs are placed at the most recent time in their placed interval. This option may be useful for plotting. The sites argument becomes arbitrary if nonstoch.bin is TRUE.
}
\value{
The output of these functions is a time-scaled tree or set of time-scaled trees, of either class phylo or multiphylo, depending on the argument ntrees. All trees are output with an element $root.time. This is the time of the root on the tree and is important for comparing patterns across trees.
}
\references{
Bapst, in prep. Time-scaling Trees of Fossil Taxa. To be submitted to Paleobiology.

Brusatte, S. L., M. J. Benton, M. Ruta, and G. T. Lloyd. 2008. Superiority, Competition, and Opportunism in the Evolutionary Radiation of Dinosaurs. Science 321(5895):1485-1488.

Hunt, G., and M. T. Carrano. 2010. Models and methods for analyzing phenotypic evolution in lineages and clades. In J. Alroy, and G. Hunt, eds. Short Course on Quantitative Methods in Paleobiology. Paleontological Society.

Lloyd, G. T., S. C. Wang, and S. L. Brusatte. 2012. Identifying Heterogeneity in Rates of Morphological Evolutio: Discrete Character Change in the Evolution of Lungfish(Sarcopterygii, Dipnoi). Evolution 66(2):330-348.

Smith, A. B. 1994. Systematics and the fossil record: documenting evolutionary patterns. Blackwell Scientific, Oxford.
}
\author{David W. Bapst, heavily inspired by code supplied by Graeme Lloyd and Gene Hunt.}

\seealso{
\code{\link{srcTimePaleoPhy}},\code{\link{binTimeData}},\code{\link{multi2di}}
}
\examples{
dev.new()
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy using the continuous range data
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)

#that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#let's add those terminal ranges back on with add.term
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#plot diversity curve 
phyloDiv(ttree)

#that tree didn't look very resolved, does it? (The curse of sampled ancestry!)
#can randomly resolve trees using the argument randres
#each resulting tree will have polytomies randomly resolved in different ways using multi2di()
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=1,randres=TRUE,add.term=TRUE,plot=TRUE)
#notice well the warning it prints!
#now let's plot the first tree (both trees will be identical because we used set.seed)
phyloDiv(ttree)
#we would need to set ntrees to a large number to get a fair sample of trees


#if we set ntrees>1, timePaleoPhy() will make multiple time-trees
ttrees<-timePaleoPhy(cladogram,rangesCont,type="basic",ntrees=9,randres=TRUE,add.term=TRUE,plot=TRUE)
#let's compare nine of them at once in a plot
layout(matrix(1:9,3,3))
for(i in 1:9){plot(ladderize(ttrees[[i]]),show.tip.label=FALSE,no.margin=TRUE)}
#they are all a bit different!
#can plot the median diversity curve with multiDiv
dev.new()
multiDiv(ttrees)

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

#using node.mins
dev.new()
#let's say we have (molecular??) evidence that node #5 is at least 1200 time-units ago
nodeDates<-rep(NA,(Nnode(cladogram)-1))
nodeDates[5]<-1200
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",
	randres=FALSE,node.mins=nodeDates,plot=TRUE)
ttree<-timePaleoPhy(cladogram,rangesCont,type="basic",
	randres=TRUE,node.mins=nodeDates,plot=TRUE)

#Using bin_timePaleoPhy to timescale with discrete interval data
#first let's use binTimeData() to bin in intervals of 1 time unit
rangesDisc<-binTimeData(rangesCont,int.length=1)
ttree<-bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,randres=TRUE,add.term=TRUE,plot=FALSE)
#notice the warning it prints!
phyloDiv(ttree)
#can also force the appearance timings not to be chosen stochastically
ttree1<-bin_timePaleoPhy(cladogram,rangesDisc,type="basic",ntrees=1,nonstoch.bin=TRUE,randres=TRUE,add.term=TRUE,plot=FALSE)
phyloDiv(ttree1)
}
