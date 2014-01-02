\name{degradeTree}
\alias{degradeTree}
\title{Randomly Collapse a Portion of Nodes on a Phylogeny}
\usage{
degradeTree(tree, prop_collapse, nCollapse = NULL, node.depth = NA,
  leave.zlb = FALSE)
}
\arguments{
  \item{tree}{A phylogeny of class 'phylo'}

  \item{prop_collapse}{Proportion of nodes to collapse}

  \item{nCollapse}{Number of nodes to collapse, can be
  supplied as an alternative to prop_collapse}

  \item{node.depth}{A number between 0 to 1, which
  conditions the depth of nodes removed. If NA, no
  conditioning (this is the default).}

  \item{leave.zlb}{If FALSE, the default option, the
  original branch length distribution is destroyed and
  branches set to zero by this function will return
  polytomies. If TRUE, then the original edge lengths are
  kept for unmodified edges, and modified edges are changed
  to zero length, and are not collapsed into polytomies.
  The removed branch length is not shifted to other edges.}
}
\value{
Returns the modified tree as an object of class phylo, with
no edge lengths by default.
}
\description{
This function removes a proportion of the total nodes in a
tree, chosen randomly, collapsing the nodes to produce a
less-resolved tree.
}
\details{
The nodes are removed at random using the basic function
sample. degradeTree can be conditioned to remove nodes of a
particular depth with greater probability/frequency by
setting node.depth to a value between zero (favoring the
removal of deep nodes close to the root) or one (shallow
nodes far from the root). Depth is evaluated based on the
number of descendant tips. If node.depth is not NA, the
relative proportion of descendants from each node is
calculated, summed to 1 and the node.depth value subtracted
from this proportion. These values are then squared,
normalized again to equal to 1 and then used as the
probabilities for sampling nodes for removal.

By default, branch lengths are removed from the input tree
prior to degradation and entirely absent from the output
tree. This is changed if argument leave.zlb is TRUE.
}
\examples{
set.seed(444)
tree <- rtree(100)
tree1 <- degradeTree(tree,0.5)

#let's compare the input and output
layout(matrix(1:2,,2))
plot(tree,show.tip.label=FALSE,use.edge.length=FALSE)
plot(tree1,show.tip.label=FALSE,use.edge.length=FALSE)

layout(1)
}
\seealso{
\code{\link{di2multi}},\code{\link{timeLadderTree}},
}

