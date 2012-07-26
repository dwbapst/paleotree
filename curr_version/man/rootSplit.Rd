\name{rootSplit}
\alias{rootSplit}
\title{Split Tip Taxa by Root Divergence}
\description{
Sorts terminal taxa into groups descended from each lineage splitting off of the root node.
}
\usage{
rootSplit(tree)
}
\arguments{
  \item{tree}{A phylo object}
}
\details{
This function can be useful for studying the descended from the first bifurcation.
}
\value{
Returns a list with each element a character vector containing the names of terminal taxa descended from each lineage splitting off of the root node.
}
\author{
David W. Bapst
}
\examples{
tree<-rtree(100)
rootSplit(tree)
}