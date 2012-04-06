\name{taxa2cladogram}
\alias{taxa2cladogram}
\title{Converting taxon data into cladogram}
\description{Convert ancestor-descendant relationships of taxa into an ideal unscaled cladogram}
\usage{taxa2cladogram(taxad, plot = F)}
\arguments{
  \item{taxad}{A five-column matrix of taxonomic data, as output by simFossilTaxa}
  \item{plot}{Should the result be plotted?}
}
\details{
This function takes a matrix that has information on the ancestor-descendant relationships of a bunch of taxa and constructs an unscaled cladogram of the hierarcially-nesting relationships among those taxa. 

The result will probably not be fully resolved, as including both ancestor and descendant taxa will make it nearly impossible to produce a fully nesting system of relationships. For example, consider a set of three morphologically-static taxa where the first is an ancestor (either direct or indirect, ala Foote, 1996) of the second and third. This set of taxa cannot be broken up into bifurcating nested relationships and will result in a polytomy for the output of this function. Any set of ancestor-descendant relationships will have many of these, as some ancestors must have more than one descendant for the clade to diversify at all.
}
\value{
The resulting phylogeny without branch lengths is output as an object of class phylo
}
\references{
Foote, M. 1996. On the Probability of Ancestors in the Fossil Record. Paleobiology 22(2):141-151.
}
\author{David W. Bapst}
\seealso{
\code{\link{simFossilTaxa}},\code{\link{taxa2phylo}},
}
\examples{
set.seed(444)
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,nExtant=0)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
cladogram<-taxa2cladogram(taxa,plot=TRUE)
plot(cladogram)
}
\keyword{datagen}