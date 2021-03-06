% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa2cladogram.R
\name{taxa2cladogram}
\alias{taxa2cladogram}
\title{Convert Simulated Taxon Data into a Cladogram}
\usage{
taxa2cladogram(taxaData, drop.cryptic = FALSE, plot = FALSE)
}
\arguments{
\item{taxaData}{A five-column matrix of taxonomic data,
as output by \code{\link{fossilRecord2fossilTaxa}} via simulations
produced using \code{\link{simFossilRecord}}. Previously, this was the
default output of the deprecated function \code{simFossilTaxa}.}

\item{drop.cryptic}{Should cryptic species be dropped (except for the
first; effectively merging the cryptic species complexes into a single apparent species)? 
\code{drop.cryptic = FALSE} by default, so cryptic species are not dropped by default.}

\item{plot}{If \code{TRUE} result the output with \code{ape::plot.phylo}.}
}
\value{
The resulting phylogeny without branch lengths is output as an
object of class \code{phylo}.

The tip labels are the rownames from the simulation input; see documentation
for \code{simFossilRecord} and \code{fossilRecord2fossilTaxa} documentation for details.
}
\description{
Convert ancestor-descendant relationships of taxa into an 'ideal' unscaled
cladogram, where taxa that could share true synapomorphies are arranged into
nested clades.
}
\details{
This function simulates an ideal cladistic process, where the relationships
of a set of morphologically static taxa is resolved into a set of nested
hierarchical relationships (a standard cladogram), as much as would be
expected given the input relationships among those taxa. taxa2cladogram uses
information on the ancestor-descendant relationships of a bunch of taxa and
constructs an unscaled cladogram of the hierarchically-nesting relationships
among those taxa. There's no actual cladistics going on, this is just a
simulation of that process. If there is any chance that a set of taxa could
be resolved into a set of nested relationships given their
ancestor-descendant relationships, they will be resolved so in the output of
taxa2cladogram. No morphological characters are considered, we just assume
that if there is a nesting relationship, then it could be resolved as such.
This makes it the "ideal" cladogram of a simulated clade.

The result will probably not be fully resolved, as including both ancestor
and descendant taxa will generally make it impossible to produce a fully
nesting system of relationships. For example, consider a set of three
morphologically-static taxa where the first is an ancestor (either direct or
indirect, ala Foote, 1996) of both the second and third. If we imagine an
ideal cladistic analysis of the morphological characters of those three
taxa, this set of taxa will be unable to be broken up into
bifurcating-nested relationships and thus result in a polytomy. Any set of
ancestor-descendant relationships will have many of these, as some ancestors
must have more than one descendant for the clade to diversify, as noted by
Wagner and Erwin, 1995.

If there are cryptic taxa present in the output from \code{simFossilRecord}, these
and any of their morphologically distinguishable descendants are collapsed
into a polytomy to simulate the expected pattern of lack of phylogenetic
resolution. In addition to this merging, cryptic taxa can be dropped via the
argument drop.cryptic, such that only the first 'species' of each cryptic
taxon assemblage is listed among the tip taxa (what we would actually expect
to obtain, as we would not recognize cryptic taxa to be treated as different OTUs). By
default, cryptic taxa are not dropped so that the same number of taxa as in
the simulated data is retained.
}
\examples{
\donttest{
set.seed(444)
record <- simFossilRecord(p = 0.1, q = 0.1, nruns = 1,
	nTotalTaxa = c(30,40), nExtant = 0)
taxa <- fossilRecord2fossilTaxa(record)
#let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
layout(1:2)
cladogram <- taxa2cladogram(taxa,plot = TRUE)
#compare the "real" time-scaled tree of taxon last occurrences (taxa2phylo) 
     #to the 'ideal' cladogram
tree <- taxa2phylo(taxa,plot = TRUE)

#testing with cryptic speciation
recordCrypt <- simFossilRecord(p = 0.1, q = 0.1, prop.cryptic = 0.5, nruns = 1,
	nTotalTaxa = c(30,40), nExtant = 0)
taxaCrypt <- fossilRecord2fossilTaxa(recordCrypt)
layout(1:2)
parOrig <- par(no.readonly = TRUE)
par(mar = c(0,0,0,0))
cladoCrypt1 <- taxa2cladogram(taxaCrypt,drop.cryptic = FALSE)
plot(cladoCrypt1)
cladoCrypt2 <- taxa2cladogram(taxaCrypt,drop.cryptic = TRUE)
plot(cladoCrypt2)

#reset plotting
par(parOrig)
layout(1) 
}
}
\references{
Foote, M. 1996 On the Probability of Ancestors in the Fossil
Record. \emph{Paleobiology} \bold{22}(2):141-151.

Wagner, P., and D. Erwin. 1995 Phylogenetic patterns as tests of speciation
models. New approaches to speciation in the fossil record. Columbia
University Press, New York:87-122.
}
\seealso{
\code{\link{simFossilRecord}}, \code{\link{taxa2phylo}}, \link{fossilRecord2fossilTaxa}
}
\author{
David W. Bapst
}
