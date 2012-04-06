\name{dropZLB}
\alias{dropZLB}
\title{Drop Terminal Zero-Length Branches}
\description{
Drops all terminal branches and tip nodes of a tree that is zero-length. Useful for paleo-trees, as ancestral taxa will often come out as these "ZLBS".
}
\usage{dropZLB(tree)}
\arguments{
  \item{tree}{A phylogeny as a phylo object}
}
\details{
This function drops tip-taxa that are attached to the tree as zero-length branches ("ZLBs"). This is advised for analyses of diversification/diversity, as these will appear as instanteous speciation+extinction events, when in reality they probably represent 'pseudo-speciation' (also referred to as anagenesis).

Note this function only drops terminal tips that are on a ZLB. If you want to drop internal zero-length branches, see \code{\link{di2multi}}

If the tree given to the function as a $root.time element, as expected for most paleo-tree objects handled by this library, that root.time is adjusted if the relative time of the root divergence changes when ZLBs are dropped.
}
\value{
Gives back a modified phylogeny as a phylo object
}
\author{David W. Bapst}
\seealso{\code{\link{phyloDiv}},\code{\link{drop.tip}}}
\examples{
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's make a tree using taxa2phylo()
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
layout(matrix(1:2,,2))
plot(ladderize(tree))
plot(ladderize(dropZLB(tree)))
#compare the two trees
}
\keyword{datagen}
