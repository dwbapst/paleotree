\name{timeSliceTree}
\alias{timeSliceTree}
\title{Taking a Timeslice of a Tree}
\description{Removes the portion of a tree after a set point in time, as if the tree had been sliced at that moment.}
\usage{timeSliceTree(ttree, sliceTime, drop.extinct = F, plot = T)}
\arguments{
  \item{ttree}{A time-scaled phylogeny of class phylo}
  \item{sliceTime}{Time to 'slice' the tree at. See details below.}
  \item{drop.extinct}{If true, drops tips that go extinct before timeSlice.}
  \item{plot}{If true, plots input and output trees for comparison.}
}
\details{
The sliceTime is always calculated as on the same as scale as ttree$root.time, which the time before present that the root divergence occurred. In other words, if root.time=100, then timeSlice=80 will slice the tree 20 time units after the root.

The function assumes that ttree will generally have an element called $root.time. If $root.time is not present as an element of ttree, then it is assumed the tip furthest from the root is at time 0 (present-day) and a new $root.time is calculated (a warning will be issued in this case).

If drop.extinct=T, then extinct tips are dropped and (if present) the $root.time of ttree is adjusted. 
}
\value{
Returns the modified phylogeny as an object of class phylo
}
\author{David W. Bapst}
\seealso{\code{\link{phyloDiv}},
}
\examples{
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
taxicDivCont(taxa[,3:4])
#that's the whole diversity curve
#with timeSliceTree we could look at the lineage accumulation curve we'd get of species sampled at a point in time
tree<-taxa2phylo(taxa)
#use timeSliceTree to make tree of relationships up until time=950 
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=FALSE)
#use drop.extinct=T to only get the tree of lineages extant at time=950
tree950<-timeSliceTree(tree,sliceTime=950,plot=TRUE,drop.extinct=TRUE)
#now its an ultrametric tree with many fewer tips...
#lets plot the lineage accumulation plot on a log scale
phyloDiv(tree950,plotLogRich=TRUE)
}
\keyword{datagen}
