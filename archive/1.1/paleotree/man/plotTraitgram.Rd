\name{plotTraitgram}
\alias{plotTraitgram}
\title{Plot a Traitgram for Continuous Traits
}
\description{This function uses maximum-likelihood ancestral trait estimation to plot a 'traitgram' (Ackerly, 2009) given a tree and a set of continuous trait values.
}
\usage{
plotTraitgram(trait, tree, trait.name = "'trait'", conf.int = T, lwd = 1.5)
}
\arguments{
  \item{trait}{a vector of continuous trait values}
  \item{tree}{a phylo object}
  \item{trait.name}{The name of the trait plotted, used for the figure's main title}
  \item{conf.int}{if true, confidence intervals are plotted}
  \item{lwd}{The line width used for the figure}
}
\details{
By default, this function uses \code{\link{ace}} from the library ape to reconstruct ancestral traits and confidence intervals with maximum likelihood. If confidence intervals are turned off, a somewhat faster function \code{\link{getAncStates}} from geiger is used.

As with many functions in the paleotree library, absolute time is always decreasing, i.e. the present day is zero.
}
\value{
Returns no value, just plots the traitgram.
}
\references{
Ackerly, D. 2009. Conservatism and diversification of plant functional traits: Evolutionary rates versus phylogenetic signal. Proceedings of the National Academy of Sciences 106(Supplement 2):19699-19706.
}
\author{David W. Bapst}
\note{
One should probably never do ancestral trait estimation without looking at the confidence intervals, as these estimates tend to be very uncertain.
}
\seealso{
\code{\link{ace}}

Also see the function \code{traitgram} in the library picante.
}
\examples{
require(geiger)
set.seed(444)
tree<-rtree(10)
trait<-rTraitCont(tree)
#first, traitgram without conf intervals
plotTraitgram(trait,tree,conf.int=FALSE)

#now, with
plotTraitgram(trait,tree)
#not much confidence, eh?
}
\keyword{hplot}
