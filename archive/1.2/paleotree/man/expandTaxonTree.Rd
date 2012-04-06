\name{expandTaxonTree}
\alias{expandTaxonTree}
\title{Extrapolating Lower Taxon Phylogenies from Higher Taxon Trees}
\description{
This tree takes a tree composed of higher taxa and a vector of lower taxa belonging to the first and produces a tree composed of the lower taxa, with the higher taxa treated as unresolved polytomies. A user can also mark higher taxa as paraphyletic and collapse them so they are not monophyetic in the outcome.
}
\usage{
expandTaxonTree(taxonTree, taxaData, collapse = NULL, plot = F)
}
\arguments{
  \item{taxonTree}{A phylo object where tips represent higher taxa}
  \item{taxaData}{A character vector of higher taxa, with elements names equal to the lower taxa.}
  \item{collapse}{A character vector of non-monophyletic higher taxa to be collapsed}
  \item{plot}{If true, plots a comparison between input and output trees}
}
\details{
The result tree will probably be a rough unresolved view of the relationships among the taxa, full of soft polytomies. This is somewhat similar to the methods used in Webb and Donoghue (2005) and Friedman (2009). Obviously, any analyses one wishes to do should be done by resolving this tree with \code{\link{multi2di}} in the ape package or the various time-scaling functions found in this package (paleotree).

The taxaData vector should have one element per lower taxon that is to be added to the tree. The names of each elements should be the desired names for the tips of the resulting lower-taxon tree that is output by this function.

The collapse vector should be used if there is doubt about the monophyly of higher taxon on the input taxontree. By default, all higher taxa are treated as monophyletic if not otherwise specified.

The output tree will always lack branch lengths and thus will not be time-scaled.
}
\value{Outputs the modified tree as an object of class phylo, with the higher taxa expanded into polytomies.}
\references{
Friedman, M. 2009. Ecomorphological selectivity among marine teleost fishes during the end-Cretaceous extinction. Proceedings of the National Academy of Sciences 106(13):5218-5223.

Webb, C. O., and M. J. Donoghue. 2005. Phylomatic: tree assembly for applied phylogenetics. Molecular Ecology Notes 5(1):181-183.
}
\author{David W. Bapst}
\seealso{
\code{\link{multi2di}},\code{\link{bind.tree}},
}
\examples{
set.seed(444)
#lets make our hypothetical simulated tree of higher taxa
taxtr<-rtree(10)
taxd<-sample(taxtr$tip.label,30,replace=TRUE)	#taxa to place within higher taxa
names(taxd)<-paste(taxd,"_x",1:30,sep="")
coll<-sample(taxtr$tip.label,3)	#what to collapse?
expandTaxonTree(taxonTree=taxtr,taxaData=taxd,collapse=coll,plot=TRUE)
}