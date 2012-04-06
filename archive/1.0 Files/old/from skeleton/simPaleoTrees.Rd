\name{simPaleoTrees}
\alias{simPaleoTrees}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simPaleoTrees(p, q, r, ntrees, nExtant = NULL, modern.samp = T, drop.zlb = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{ntrees}{
%%     ~~Describe \code{ntrees} here~~
}
  \item{nExtant}{
%%     ~~Describe \code{nExtant} here~~
}
  \item{modern.samp}{
%%     ~~Describe \code{modern.samp} here~~
}
  \item{drop.zlb}{
%%     ~~Describe \code{drop.zlb} here~~
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
function (p, q, r, ntrees, nExtant = NULL, modern.samp = T, drop.zlb = T) 
{
    require(ape)
    res <- rmtree(ntrees, 2)
    for (i in 1:ntrees) {
        rerun <- T
        while (rerun) {
            taxa <- simFossilTaxa(p, q, w = 0, u = 0, nruns = 1, 
                mintaxa = 2, maxtaxa = 500, maxtime = 1000, nExtant = nExtant, 
                plot = T)
            ranges <- sampRanges(taxa, r, min.taxa = 0, modern.samp = modern.samp)
            if (sum(!is.na(ranges[, 1])) > 1) {
                tree <- taxa2phylo(taxa, obs_time = ranges[, 
                  2], plot = T)
                if (drop.zlb == T) {
                  tree <- dropZLB(tree)
                }
                if (all(!is.na(tree))) {
                  rerun <- F
                }
            }
        }
        res[[i]] <- tree
        if ((i\%\%(ntrees/10)) == 0) {
            cat(round(i * 100/ntrees), "\% ", sep = "")
        }
    }
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
