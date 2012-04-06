\name{dropZLB}
\alias{dropZLB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dropZLB(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
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
function (tree) 
{
    require(ape)
    dropZLB_onetree <- function(tree) {
        drop_e <- (tree$edge[, 2] < (Ntip(tree) + 1)) & (tree$edge.length == 
            0)
        drop_t <- (tree$edge[, 2])[drop_e]
        if ((Ntip(tree) - length(drop_t)) > 1) {
            tree1 <- drop.tip(tree, drop_t)
            orig_dist <- dist.nodes(tree)[which(tree1$tip.label[1] == 
                tree$tip.label), Ntip(tree) + 1]
            new_dist <- dist.nodes(tree1)[1, Ntip(tree1) + 1]
            tree1$root.time <- tree$root.time - (orig_dist - 
                new_dist)
            res <- tree1
        }
        else {
            res <- NA
        }
        res
    }
    if (any(names(tree) == "edge")) {
        res <- dropZLB_onetree(tree)
    }
    else {
        res <- lapply(tree, function(x) dropZLB_onetree(x))
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
