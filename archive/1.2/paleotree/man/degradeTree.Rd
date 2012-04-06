\name{degradeTree}
\alias{degradeTree}
\title{Randomly collapse nodes on a phylogeny}
\description{This function removes a proportion of the total nodes in a tree, chosen randomly, collapsing the nodes such that it causes the tree to be more unresolved.}
\usage{degradeTree(tree, prop_collapse, node.depth = NA)}
\arguments{
  \item{tree}{A tree of class phylo}
  \item{prop_collapse}{Proportion of nodes to collapse}
  \item{node.depth}{A number between 0 to 1, which conditions the depth of nodes removed. If NA, no conditioning (this is the default).}
}
\details{
The nodes are removed at random using sample(). This function can be conditioned to remove nodes of a particular depth with greater probability/frequency by setting node.depth to some value between zero (deep nodes close to the root) or one (shallow nodes far from the root). Depth is evaluated based on the number of descendant tips. If node.depth is not NA, the relative proportion of descendants from each node is calculated, summed to 1 and the node.depth value subtracted from this proportion. These values are then squared, normalized again to equal to 1 and then used as the probabilities for sampling nodes for removal.

Branch lengths are removed from the input tree prior to degradation and entirely absent from the output tree.
}
\value{
Returns the modified tree as an object of class phylo, with no edge lengths.
}
\author{David W. Bapst}

\seealso{\code{\link{di2multi}}}
\examples{
set.seed(444)
tree<-rtree(100)
tree1<-degradeTree(tree,0.5)
#let's compare the input and output
layout(matrix(1:2,,2))
plot(tree,show.tip.label=FALSE,use.edge.length=FALSE)
plot(tree1,show.tip.label=FALSE,use.edge.length=FALSE)
}