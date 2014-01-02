\name{timeLadderTree}
\alias{timeLadderTree}
\title{Resolve Polytomies by Order of First Appearance}
\usage{
timeLadderTree(tree, timeData)
}
\arguments{
  \item{tree}{A phylo object}

  \item{timeData}{Two-column matrix of per-taxon first and
  last occurrances in absolute continous time}
}
\value{
Returns the modified tree as an object of class phylo, with
no edge lengths.
}
\description{
Resolves polytomies in trees with lineages arranged in a
pectinate pattern (i.e. a ladder-like subtree), ordered by
the time of first appearance (FAD) for each lineage.
}
\details{
This method of resolving polytomies assumes that the order
of stratigraphic appearance perfectly depicts the order of
branching. This may not be a good assumption for poorly
sampled fossil records.

This function is for resolving trees when a continuous
time-scale is known. For discrete time-scales, see the
function bin_timePaleoPhy.

Taxa with the same identical first appearance date will be
ordered randomly. Thus, the output is slightly stochastic,
but only when ties exist. This is probably uncommon with
real data on continuous time-scales.

Taxa not shared between the tree and the timeData matrix,
or listed as having a FAD or LAD of NA in timeData will be
dropped and will not be included in the output tree.

See this blog post for more information:
http://nemagraptus.blogspot.com/2012/07/resolving-polytomies-according-to.html
}
\examples{
set.seed(444)
taxa<-simFossilTaxa(0.1,0.1,mintaxa=100)
tree<-taxa2cladogram(taxa)
ranges<-sampleRanges(taxa,r=0.5)
tree1<-timeLadderTree(tree,ranges)
layout(1:2)
plot(ladderize(tree),show.tip.label=FALSE)
plot(ladderize(tree1),show.tip.label=FALSE)

#an example with applying timeLadderTree to discrete time data
rangeData<-binTimeData(ranges,int.len=5)	#sim discrete range data
tree2<-bin_timePaleoPhy(tree,timeList=rangeData,timeres=TRUE)
plot(ladderize(tree),show.tip.label=FALSE)
plot(ladderize(tree2),show.tip.label=FALSE)

layout(1)
}
\author{
David W. Bapst
}
\seealso{
\code{\link{di2multi}}
}

