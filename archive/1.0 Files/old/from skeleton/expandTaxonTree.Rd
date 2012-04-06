\name{expandTaxonTree}
\alias{expandTaxonTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
expandTaxonTree(taxonTree, taxaData, collapse = NULL, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxonTree}{
%%     ~~Describe \code{taxonTree} here~~
}
  \item{taxaData}{
%%     ~~Describe \code{taxaData} here~~
}
  \item{collapse}{
%%     ~~Describe \code{collapse} here~~
}
  \item{plot}{
%%     ~~Describe \code{plot} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (taxonTree, taxaData, collapse = NULL, plot = T) 
{
    require(ape)
    tree <- taxonTree
    tree$edge.length <- rep(1, Nedge(tree))
    for (i in unique(taxaData)) {
        tip <- which(tree$tip.label == i)
        if (length(collapse) > 0) {
            if (any(collapse == i)) {
                tree$edge.length[which(tree$edge[, 2] == tip)] <- 0
            }
        }
        cotaxa <- names(taxaData)[taxaData == i]
        repTree <- stree(length(cotaxa))
        repTree$edge.length <- rep(1, length(cotaxa))
        repTree$tip.label <- cotaxa
        tree <- bind.tree(tree, repTree, tip)
    }
    tree1 <- di2multi(tree)
    tree1$edge.length <- NULL
    tree1 <- collapse.singles(tree1)
    tree1 <- read.tree(text = write.tree(tree1))
    if (plot == T) {
        layout(matrix(1:2, , 2))
        plot(taxonTree)
        plot(tree1)
    }
    return(tree1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
