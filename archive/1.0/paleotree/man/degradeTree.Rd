\name{degradeTree}
\alias{degradeTree}
\title{Randomly collapse nodes on a phylogeny}
\description{This function removes a proportion of the total nodes in a tree, chosen randomly, collapsing the nodes such that it causes the tree to be more unresolved.}
\usage{degradeTree(tree, prop_collapse, node.depth = NA)}
\arguments{
  \item{tree}{A tree of class phylo}
  \item{prop_collapse}{Proportion of nodes to collapse}
  \item{node.depth}{Conditions depth of nodes removed. Default is NA (no conditioning).}
}
\details{
This function can be conditioned to remove nodes of a particular depth with greater probability by setting node.depth to some value between zero (deep nodes close to the root) or one (shallow nodes far from the root). 
}
\value{
Returns the modified tree as an object of class phylo
}
\author{David W. Bapst}

\seealso{\code{\link{di2multi}}}
\examples{
set.seed(444)
tree<-rtree(100)
tree1<-degradeTree(tree,0.5)
layout(matrix(1:2,,2))
plot(tree,show.tip.label=FALSE,use.edge.length=FALSE)
plot(tree1,show.tip.label=FALSE,use.edge.length=FALSE)
}
\keyword{datagen}