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
Returns the modified phylogeny as an object of class phylo
}
\author{David W. Bapst}
\seealso{
\code{\link{speciationalTree}},\code{\link{timePaleoPhy}},
}
\examples{
set.seed(444)
tree<-rtree(10)
layout(matrix(1:2,,2))
plot(tree);plot(unitLengthTree(tree))
}
\keyword{datagen}

