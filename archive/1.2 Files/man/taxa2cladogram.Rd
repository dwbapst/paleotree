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
This function simulates an ideal cladistic process, where the relationships of a set of morphologically static taxa is resolved into a set of nested hierarchial relationships (a standard cladogram), as much as would be expected given the input relationships among those taxa. taxa2cladogram uses information on the ancestor-descendant relationships of a bunch of taxa and constructs an unscaled cladogram of the hierarcially-nesting relationships among those taxa. There's no actual cladistics going on, this is just a simulation of that process. If there is any chance that a set of taxa could be resolved into a set of nested relationships given their ancestor-descendant relationships, they will be resolved so in the output of taxa2cladogram. No morphological characters are considered, we just assume that if there is a nesting relationship, then it could be resolved as such. This makes it the "ideal" cladogram of a simulated clade.

The result will probably not be fully resolved, as including both ancestor and descendant taxa will generally make it impossible to produce a fully nesting system of relationships. For example, consider a set of three morphologically-static taxa where the first is an ancestor (either direct or indirect, ala Foote, 1996) of both the second and third. If we imagine an ideal cladistic analysis of the morphological characters of those three taxa, this set of taxa will be unable to be broken up into bifurcating-nested relationships and thus result in a polytomy. Any set of ancestor-descendant relationships will have many of these, as some ancestors must have more than one descendant for the clade to diversify.
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
taxa<-simFossilTaxa(p=0.1,q=0.1,nruns=1,mintaxa=20,maxtaxa=30,maxtime=1000,maxExtant=0)
#let's use taxa2cladogram() to get the 'ideal' cladogram of the taxa
layout(matrix(1:2,,2))
cladogram<-taxa2cladogram(taxa,plot=TRUE)
#compare the real relationships (taxa2phylo) to the ideal cladogram
tree<-taxa2phylo(taxa,plot=TRUE)
}