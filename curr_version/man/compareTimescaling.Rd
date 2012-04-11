\name{compareTimescaling}
\alias{compareTimescaling}
\alias{compareNodeAges}
\alias{compareTermBranches}
\title{Comparing the Time-Scaling of Trees}
\description{These functions take two trees and then calculate the changes in node ages (for compareNodeAges) for shared clades or terminal branch lengths leading to shared tip taxa (for compareTermBranches).}
\usage{
compareNodeAges(tree1,tree2)

compareTermBranches(tree1,tree2)
}
\arguments{
  \item{tree1}{A time-scaled phylogeny of class 'phylo'}
  \item{tree2}{A time-scaled phylogeny of class 'phylo'}
}
\details{
These functions are for comparing the time-scaling of two trees. Any taxa not-shared on both trees are dropped before analysis, based on tip labels. 

As with many paleotree functions, calculations relating to are done with respect to any included $root.time elements. If these are not present, the latest tip is assumed to be at the present day (time=0).

compareNodeAges calculates the changes in the clade ages among those clades shared by the two trees, relative to the first tree in absolute time. For example, a shift of +5 means the clade originates 5 time units later in absolute on the second tree, while a shift of -5 means the clade originated 5 time-units prior on the second tree.

compareTermBranches calculates the changes in the terminal branch lengths attached to tip taxa shared by the two trees, relative to the first tree. Thus, a shift of +5 means that terminal taxon has a terminal branch which is five time-units longer.
}
\value{
Returns a vector of temporal shifts for nodes or terminal branches, dependent on the function applied. If compareTermBranches is used, the vector will have the tip names as labels.
}
\author{David W. Bapst}
\seealso{
\code{\link{taxa2phylo}},\code{\link{phyloDiv}}
}
\examples{
dev.new()
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#get the true tree
tree1<-taxa2phylo(taxa)
#simulate a fossil record with imperfect sampling with sampleRanges()
rangesCont<-sampleRanges(taxa,r=0.5)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#Now let's try timePaleoPhy using the continuous range data
tree2<-timePaleoPhy(cladogram,rangesCont,type="basic")
#let's look at the distribution of node shifts
hist(compareNodeAges(tree1,tree2))
#let's look at the distribution of terminal branch lengths
hist(compareTermBranches(tree1,tree2))

set.seed(444)
tree1<-rtree(10)
tree2<-tree1
tree1$root.time<-10
compareNodeAges(tree1,tree2)
compareTermBranches(tree1,tree2)
}