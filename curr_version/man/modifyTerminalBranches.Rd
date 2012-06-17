\name{modifyTerminalBranches}
\alias{dropZLB}
\alias{dropExtinct}
\alias{dropExtant}
\alias{addTermBranchLength}
\title{Modify or Drop Terminal Branches of Various Types}
\description{
These functions modify terminal branches or drop certain terminal branches based on various criteria.
}
\usage{
dropZLB(tree)
dropExtant(tree, tol = 0.01)
dropExtinct(tree, tol = 0.01, ignore.root.time = FALSE)
addTermBranchLength(tree, addtime = 0.001)
}
\arguments{
  \item{tree}{A phylogeny as a phylo object}
  \item{tol}{Tolerance for determining modern age; used for distinguishing extinct from extant taxa. Tips which end within 'tol' of the furthest distance from the root will be treated as 'extant' taxa for the purpose of keeping or dropping.}
  \item{ignore.root.time}{Ignore root.time in calculating which tips are extinct? root.time will still be adjusted}
  \item{addtime}{Extra amount of time to add to all terminal branch lengths.}
}
\details{
DropZLB drops tip-taxa that are attached to the tree via zero-length terminal branches ("ZLBs"). This is sometimes useful for paleo-trees, as various time-scaling methods often produce these ZLBs, taxa whose early appearance causes them to be functionally interpreted as ancestors in some time-scaling methods. Removing ZLBs is advised for analyses of diversification/diversity, as these will appear as simultaneous speciation/extinction events. Note this function only drops tips attached to a terminal zero-length branch; if you want to collapse internal zero-length branches, see the ape function \code{\link{di2multi}}.

DropExtinct drops all terminal branches which end before the modern (i.e. extinct taxa). DropExtant drops all terminal branches which end at the modern (i.e. extant/still-living taxa). In both cases, the modern is defined based on tree$root.time if available, or the modern is inferred to be the point in time when the tip furthest from the root (the latest tip) terminates.

If the input tree has a $root.time element, as expected for most paleo-tree objects handled by this library, that root.time is adjusted if the relative time of the root divergence changes when terminal branches are dropped. Adjusted root.times are only given if the input tree has root.times.

addTermBranchLength adds an amount equal to the argument 'addtime' to the terminal branch lengths of the tree. If there is a $root.time element, this is increased by an amount equal to addtime. A negative amount can be input to reduce the length of terminal branches. However, if negative branch lengths are produced, the function fails and a warning is produced.
}
\value{
Gives back a modified phylogeny as a phylo object
}
\author{David W. Bapst}
\seealso{\code{\link{phyloDiv}}, \code{\link{drop.tip}}, \code{\link{compareTermBranches}}}
\examples{
set.seed(444)
#Simulate some fossil ranges with simFossilTaxa
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#simulate a fossil record with imperfect sampling with sampleRanges
rangesCont<-sampleRanges(taxa,r=0.5)
#Now let's make a tree using taxa2phylo
tree <- taxa2phylo(taxa,obs_time=rangesCont[,2])
#compare the two trees
layout(1:2)
plot(ladderize(tree))
plot(ladderize(dropZLB(tree)))
layout(1)

#example using dropExtinct and dropExtant
set.seed(444)
taxa <- simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=40,maxtime=1000,maxExtant=20)
tree<-taxa2phylo(taxa)
phyloDiv(tree)
tree1 <- dropExtinct(tree)
phyloDiv(tree1)
tree2 <- dropExtant(tree)
phyloDiv(tree2)

#example using addTermBranchLength
set.seed(444)
treeA <- rtree(10)
treeB <- addTermBranchLength(treeA,1)
compareTermBranches(treeA,treeB)
}
