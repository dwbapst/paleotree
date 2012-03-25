\name{timeSliceTree}
\alias{timeSliceTree}
\title{Taking a Timeslice of a Tree}
\description{Removes the portion of a tree after a set point in time, as if the tree after that moment had been sliced away.}
\usage{timeSliceTree(ttree, sliceTime, drop.extinct = F, plot = T)}
\arguments{
  \item{ttree}{A time-scaled phylogeny of class phylo}
  \item{sliceTime}{Time to 'slice' the tree at. See details below.}
  \item{drop.extinct}{If true, drops tips that go extinct before timeSlice.}
  \item{plot}{If true, plots input and output trees for comparison.}
}
\details{
The function assumes that ttree will generally have an element called $root.time, which is the time before present that the root divergence occurred. If $root.time is not present as an element of ttree, then it is assumed the tip furthest from the root is at time 0 (present-day) and a new $root.time is calculated (a warning will be issued in this case).

The sliceTime is always calculated as on the same scale as ttree$root.time. In other words, if root.time=100, then timeSlice=80 will slice the tree 20 time units after the root.

If drop.extinct=T, then extinct tips are dropped and (if present) the $root.time of ttree is adjusted. This is done using the function dropExtinct. 
}
\value{
Returns the modified phylogeny as an object of class phylo. Tip labels for cut branches which held multiple descendant tips will be the label for the earliest appearing tip descendant of that branch.
}
\author{David W. Bapst}
\seealso{\code{\link{phyloDiv}},\code{\link{dropExtinct}},\code{\link{dropExtant}}

Also see the function treeSlice in the library phytools, which will slice a tree at some point in and return all the subtrees which remain after the slicing time. (Effectively the reverse of timeSliceTree.)
}
\examples{
dev.new()
#a neat example of using phyloDiv with timeSliceTree 
 #to simulate doing molecular-phylogeny studies 
 #of diverification...in the past
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
taxicDivCont(taxa)
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve 
 #we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=TRUE)
}