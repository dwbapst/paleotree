\name{unitLengthTree}
\alias{unitLengthTree}
\title{Scale Tree to Unit-Length}
\description{Rescales all edges of a phylogeny to be equal to 1 ("unit-length").}
\usage{unitLengthTree(tree)}
\arguments{
  \item{tree}{an object of class phylo}
}
\details{
No additional details.
}
\value{
Returns the modified phylogeny as an object of class phylo. Any $root.time element is removed.
}
\author{David W. Bapst}
\seealso{
As an alternative to using unitLengthTree in comparative studies, see \code{\link{timePaleoPhy}}

See also \code{speciationalTree} in the package geiger, which does essentially the same thing as unitLengthTree.
}
\examples{
set.seed(444)
tree <- rtree(10)
layout(1:2)
plot(tree)
plot(unitLengthTree(tree))
layout(1)
}
