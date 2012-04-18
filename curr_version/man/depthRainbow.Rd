\name{depthRainbow}
\alias{depthRainbow}
\title{Paint Tree Branch Depth by Color}
\description{
Paints the edges of a phylogeny with colors relative to their depth.
}
\usage{
depthRainbow(tree)
}
\arguments{
  \item{tree}{A phylo object}
}
\details{
The only purpose of this function is to make an aesthetically-pleasing graphic of one's tree, where branches are color-coded, using a rainbow palette, relative to their depth. Depth is defined relative to the number of branching nodes between the basal node of a branch and the root, not the absolute distance (i.e. branch length) to the root or the distance from the tips.
}
\value{
No value returned, just plots a colorful phylogeny.
}
\author{
David W. Bapst
}
\examples{
set.seed(444)
tree <- rtree(500)
depthRainbow(tree)
}