\name{simFossilTaxa_SRCond}
\alias{simFossilTaxa_SRCond}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
simFossilTaxa_SRCond(r, p, q, w = 0, u = 0, nruns = 1, avgtaxa, maxtime = 100, nExtant = NULL, plot = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{u}{
%%     ~~Describe \code{u} here~~
}
  \item{nruns}{
%%     ~~Describe \code{nruns} here~~
}
  \item{avgtaxa}{
%%     ~~Describe \code{avgtaxa} here~~
}
  \item{maxtime}{
%%     ~~Describe \code{maxtime} here~~
}
  \item{nExtant}{
%%     ~~Describe \code{nExtant} here~~
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
function (r, p, q, w = 0, u = 0, nruns = 1, avgtaxa, maxtime = 100, 
    nExtant = NULL, plot = T) 
{
    N <- avgtaxa/(1 - exp(-r/q))
    results <- simFossilTaxa(p, q, w, u, nruns, mintaxa = (N), 
        maxtaxa = (N + (N/2)), maxtime, nExtant, plot)
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
