\name{droppingBranches}
\alias{dropZLB}
\alias{dropExtinct}
\alias{dropExtant}
\title{Drop Terminal Branches of Various Types}
\description{
Drops all terminal branches and tip nodes of a tree that is zero-length. Useful for paleo-trees, as ancestral taxa will often come out as these "ZLBS".
}
\usage{
dropZLB(tree)
dropExtant(tree, tol = 0.01)
dropExtinct(tree, tol = 0.01, ignore.root.time = F)
}
\arguments{
  \item{tree}{A phylogeny as a phylo object}
  \item{tol}{Tolerance of determining modern age for distinguishing extinct from extant taxa}
  \item{ignore.root.time}{Ignore root.time in calculating which tips are extinct? root.time will still be adjusted}
}
\details{
DropZLB drops tip-taxa that are attached to the tree as zero-length branches ("ZLBs"). This is advised for analyses of diversification/diversity, as these will appear as instanteous speciation+extinction events, when in reality they probably represent 'pseudo-speciation' (also referred to as anagenesis).Note this function only drops terminal tips that are on a ZLB. If you want to drop internal zero-length branches, see \code{\link{di2multi}}

DropExtinct drops all terminal branches which end before the modern (i.e. extinct taxa). DropExtant drops all terminal branches which end at the modern (i.e. extant/still-living taxa). In both cases, the modern is defined based on tree$root.time if available, or the modern is inferred to when the tip furthest from the root terminates in tree.

If the tree given to the function as a $root.time element, as expected for most paleo-tree objects handled by this library, that root.time is adjusted if the relative time of the root divergence changes when terminal branches are dropped. Adjusted root.times are only given if the input tree has root.times.
}
\value{
Gives back a modified phylogeny as a phylo object
}
\author{David W. Bapst}
\seealso{\code{\link{phyloDiv}},\code{\link{drop.tip}}}
\examples{
##Simulate some fossil ranges with simFossilTaxa()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's make a tree using taxa2phylo()
tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
#compare the two trees
layout(matrix(1:2,,2))
plot(ladderize(tree))
plot(ladderize(dropZLB(tree)))


#example using dropExtinct and dropExtant
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=40,maxtime=1000,maxExtant=20)
tree<-taxa2phylo(taxa)
phyloDiv(tree)
tree1<-dropExtinct(tree)
phyloDiv(tree1)
tree2<-dropExtant(tree)
phyloDiv(tree2)
graphics.off()
}
