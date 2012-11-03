\name{candleTaxa}
\alias{candleTaxa}
\alias{simCandleTaxa}
\alias{trueCandle}
\alias{deadTree}
\title{Simulating Extinct Clades of Monophyletic Taxa}
\description{
This function simulates the diversification of clades composed of monophyletic taxa, which are distinguished in a fashion completely alternative to way taxa are defined in the simulation functions simFossilTaxa, taxa2cladogram and taxa2phylo.
}
\usage{
deadTree(ntaxa, sumRate=0.2)
simCandleTaxa(ntaxa, sumRate=0.2)
trueCandle(candleRes, time.obs)
}
\arguments{
  \item{ntaxa}{Number of monophyletic taxa (tip terminals) to be included on the simulated tree}
  \item{sumRate}{The sum of the instantaneous branching and extinction rates; see below.}
  \item{candleRes}{The list output produced by simCandleTaxa}
  \item{time.obs}{A per-taxon vector of times of observation for the taxa in candleRes}
}
\details{
deadTree generates a time-scaled topology for an entirely extinct clade of a specific number of tip taxa. Because the clade is extinct and assumed to have gone extinct in the distant past, many details of typical birth-death simulators can be ignored. If a generated clade is already conditioned upon the (a) that some number of taxa was reached and (b) then the clade went extinct, the topology (i.e. the distribution of branching and extinction events) among the branches should be independent of the actual generating rate. The frequency of nodes is a simple mathematical function of the number of taxa (i.e. number of nodes is the number of taxa -1) and their placement should completely random, given that we generally treat birth-death processes as independent Poisson processes. Thus, in terms of generating the topology, this function is nothing but a simple wrapper for the ape function rtree, which randomly places splits among a set of taxa using a simple algorithm (see Paradis, 2012). To match the expectation of a birth-death process, new branch lengths are calculated as an exponential distribution with mean 1/sumRate, where sumRate represents the sum of the branching and extinction rates. Although as long as both the branching rate and extinction rates are more than zero, any non-ultrametric tree is possible, only when the two rates are non-zer and equal to each other will there be a high chance of getting an extinct clade with many tips. Any analyses one could do on a tree such as this will almost certainly give estimates of equal branching and extinction rates, just because all taxa are extinct.

simCandleTaxa produces 'candle' datasets; datasets of clades where the set of distinguishable taxa are defined as intrinsically monophyletic. (The 'candle' names comes from a horticultural practice which the package maintainer thinks of as being conceptually similar.) On theoretical terms, candle datasets are what would occur if (a) only descendant lineages can be sample and (b) all taxa are immediately differentiated as of the last speciation event and continue to be so differentiated until they go extinct. In practice, this means the taxa on such a tree would represent a sample of all the terminal branches, which start with some speciation event and end in an extinction event. These are taken to be the true original ranges of these taxa. No further taxa can be sampled than this set, whatsoever. Note that the differentiation here is a result of a posteriori consideration of the phylogeny: one can't even know what lineages could be sampled or the actual start points of such taxa until after the entire phylogeny of a group of organisms is generated. simCandleTaxa uses deadTree to make a phylogeny, so the only datasets produced are of extinct clades.

Because all evolutionary history prior to any branching events is unsampled, this model is somewhat agnostic about the general model of differentiation among lineages. The only thing that can be said is that synapomorphies are assumed to be potentially present along every single branch, such that in an ideal scenario every clade could be defined. This would suggest very high anagenesis or bifurcation.

Because the set of observable taxa is a limited subset of the true evolution history, the true taxon ranges are not a faithful reproduction of the true diversity curve. See an example below.

trueCandle is analagous to the function of taxa2phylo, in that it outputs the time-scaled-phylogeny for a candle-taxon dataset for some times of observations. Unlike with the use of taxa2phylo on the output on simFossilTaxa, there is no need to use trueCandle to obtain the true phylogeny when times of extinction are the times of observation; just get the $tree element from the result output by simCandleTaxa. 
Also unlike with taxa2phylo, the cladistic topology of relationships among morphotaxa never changes as a function of time of observation. For obtaining the 'ideal cladogram' of relationships among the candle taxa, merely take the $tree element of the output from simCandleTaxa and remove the branch lengths (see below for an example).

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.
}
\value{
deadTree gives time-scaled phylo object, with a $root.time element. As discussed above, the result is always an extinct phylogeny of exactly ntaxa.

simCandleTaxa produces a list with two components: $taxonRanges which is a two-column matrix where each row gives the true first and last appearance of observable taxa and $tree which is a time-scaled phylogeny with end-points at the true last appearance time of taxa.

trueCandle produces a time-scaled tree as a phylo object, which describes the relationships of populations at the times of observation given in the time.obs argument.
}
\references{
Paradis, E. (2012) \emph{Analysis of Phylogenetics and Evolution with R (Second Edition).} New York: Springer.
}
\author{
David W. Bapst
}
\seealso{
deadtree is simply a wraper of the function \code{rtree} in ape.

For a very different way of simulating diversification in the fossil record, see \code{\link{simFossilTaxa}}, \code{\link{taxa2phylo}} and \code{\link{taxa2cladogram}}.
}
\examples{
set.seed(444)
#example for 20 taxa
candleRes<-simCandleTaxa(20)

#let look at the taxa...
taxa<-candleRes$taxonRanges
taxicDivCont(taxa)
#because ancestors don't even exist as taxa
	#the true diversity curve can go to zero
	#kinda bizarre!

#the tree should give a better idea
tree<-candleRes$tree
phyloDiv(tree)
#well, okay, its a tree. 

#get the 'ideal cladogram' ala taxa2cladogram
    #much easier with candle-simulations as no paraphyletic taxa
cladogram<-tree
cladogram$edge.length<-NULL
plot(cladogram)

#trying out trueCandle
#random times of observation: uniform distribution
time.obs<-apply(taxa,1,function(x) runif(1,x[2],x[1]))
tree1<-trueCandle(candleRes,time.obs)
layout(1:2)
plot(tree)
plot(tree1)

#let's look at the change in the terminal branches
plot(tree$edge.length,tree1$edge.length)
#can see some edges are shorter on the new tree, cool

#let's now simulate sampling and use FADs
layout(1:2)
plot(tree);axisPhylo()
FADs<-sampleRanges(candleRes$taxonRanges,r=0.1)[,1]
tree1<-trueCandle(candleRes,FADs)
plot(tree1);axisPhylo()

#can condition on sampling some average number of taxa
#analagous to simFossilTaxa_SRcond
r<-0.1
avgtaxa<-50
sumRate<-0.2
ntaxa_orig<-avgtaxa/(r/(r+sumRate))	#avg number necc for an avg number sampled
candles<-simCandleTaxa(ntaxa=ntaxa_orig,sumRate=sumRate)
#note that conditioning must be conducted using full sumRate
#this is because durations are functions of both rates
#just like in bifurcation

layout(1)
}