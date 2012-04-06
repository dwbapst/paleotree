\name{degradeTree}
\alias{degradeTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
degradeTree(tree, prop_collapse, node.depth = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{prop_collapse}{
%%     ~~Describe \code{prop_collapse} here~~
}
  \item{node.depth}{
%%     ~~Describe \code{node.depth} here~~
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
function (tree, prop_collapse, node.depth = NA) 
{
    tree <- collapse.singles(tree)
    edge <- (1:length(tree$edge))[which(tree$edge[, 2] > Ntip(tree))]
    if (is.na(node.depth)) {
        cedge <- sample(edge, round(prop_collapse * length(edge)))
    }
    else {
        node_pdesc <- sapply(prop.part(tree), length)/Ntip(tree)
        edge_pdesc <- node_pdesc[tree$edge[edge, 2] - Ntip(tree)]
        edge_prob <- (edge_pdesc - node.depth)^2
        edge_prob <- edge_prob/sum(edge_prob)
        cedge <- sample(edge, round(prop_collapse * length(edge)), 
            prob = edge_prob)
    }
    tree$edge.length <- rep(1, Nedge(tree))
    tree$edge.length[cedge] <- 0
    tree <- di2multi(tree)
    tree <- collapse.singles(tree)
    tree$edge.length <- NULL
    return(tree)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
